#ifndef APM_APM_HPP_
#define APM_APM_HPP_

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <numbers>
#include <optional>
#include <queue>
#include <ranges>
#include <vector>

namespace apm::details {

struct point {
  double x;
  double y;

  constexpr point(double x, double y)
      : x(x)
      , y(y) {}

  /**
   * Returns the hypervolume of this point w.r.t. a reference point `r`.
   */
  [[nodiscard]] constexpr auto hv(point const &r) const -> double {
    return (x - r.x) * (y - r.y);
  }
};

/**
 * Returns whether two doubles are approximately equal.
 *
 * This uses the formula |a - b| <= atol + rtol * |b| taken from
 * numpy.isclose,
 * https://numpy.org/doc/stable/reference/generated/numpy.isclose.html
 */
[[nodiscard]] constexpr auto is_close(double const &a, double const &b, double rtol = 1e-8, double atol = 1e-30) -> bool {
  return std::abs(a - b) <= atol + rtol * std::abs(b);
}

[[nodiscard]] constexpr auto is_point_close(point const &a, point const &b, double rtol = 1e-8, double atol = 1e-30) -> bool {
  return is_close(a.x, b.x, rtol, atol) && is_close(a.y, b.y, rtol, atol);
}

struct segment {
  point a;
  point b;

  constexpr segment(point a, point b)
      : a(a)
      , b(b) {
    assert(a.x < b.x);
    assert(a.y > b.y);
  }

  /**
   * Returns the hypervolume of this segment w.r.t. a reference point `r`.
   */
  [[nodiscard]] constexpr auto hv(point const &r) const -> double {
    assert(a.x >= r.x);
    assert(a.y >= r.y);
    assert(b.x >= r.x);
    assert(b.y >= r.y);
    return a.hv(r) + b.hv(r) + (b.x - a.x) * (a.y - b.y) / 2.0;
  }

  /**
   * Returns the point in the segment that gives maximal hypervolume
   * w.r.t. a reference point `r`.
   *
   * Assumptions: segment dominates reference point.
   */
  [[nodiscard]] constexpr auto max_hv_point(point const &r) const -> point {
    assert(a.x >= r.x);
    assert(a.y >= r.y);
    assert(b.x >= r.x);
    assert(b.y >= r.y);

    auto const lm = (b.y - a.y) / (b.x - a.x);
    auto const lb = a.y - lm * a.x;
    auto const cx = (r.y - lb) / lm;
    auto const zx = r.x + (cx - r.x) / 2.0;
    if (zx > b.x || is_close(zx, b.x)) {
      return b;
    } else if (zx < a.x || is_close(zx, a.x)) {
      return a;
    } else {
      auto const cy = lm * r.x + lb;
      auto const zy = r.y + (cy - r.y) / 2.0;
      return point{zx, zy};
    }
  }
};

/**
 * Generate `n` piecewise approximation segments given curvature parameter `d`.
 */
[[nodiscard]] auto piecewise_segments(size_t n, double d) -> std::vector<segment> {
  auto res = std::vector<segment>{};
  res.reserve(n);

  auto const pi2 = std::numbers::pi_v<double> / 2.0;
  auto const p = 2.0 / d;
  auto const s = pi2 / static_cast<double>(n);
  auto a = point{0.0, 1.0};
  for (size_t i = 1; i < n; ++i) {
    auto const te = std::min(s * static_cast<double>(i), pi2);
    auto const b = point{std::pow(std::sin(te), p), std::pow(std::cos(te), p)};
    res.emplace_back(a, b);
    a = b;
  }
  res.emplace_back(a, point{1.0, 0.0});
  return res;
}

class exact_region {
 private:
  std::vector<segment> m_segments;
  point m_ref;
  std::optional<point> m_opt_point1;
  std::optional<point> m_opt_point2;
  double m_opt_point_hv;

 public:
  using segments_type = std::vector<segment>;

  template <typename Segments>
  constexpr exact_region(Segments &&segments, point ref)
      : m_segments(std::forward<Segments>(segments))
      , m_ref(ref)
      , m_opt_point1(std::nullopt)
      , m_opt_point2(std::nullopt)
      , m_opt_point_hv(0.0) {
    if (m_segments.empty())
      return;
    m_opt_point1 = m_segments.front().max_hv_point(m_ref);
    m_opt_point_hv = m_opt_point1->hv(m_ref);
    for (auto it = std::next(m_segments.begin()); it != m_segments.end(); ++it) {
      auto const p = it->max_hv_point(m_ref);
      auto const hv = p.hv(m_ref);
      if (is_close(hv, m_opt_point_hv)) {
        m_opt_point2 = p;
      } else if (hv > m_opt_point_hv) {
        m_opt_point_hv = hv;
        m_opt_point1 = p;
        m_opt_point2 = std::nullopt;
      }
    }
  }

  [[nodiscard]] constexpr auto operator<(exact_region const &o) const -> bool {
    return m_opt_point_hv < o.m_opt_point_hv;
  }

  [[nodiscard]] constexpr auto max_hv() const {
    return m_opt_point_hv;
  }

  [[nodiscard]] auto empty() const -> bool {
    return m_segments.empty();
  }

  [[nodiscard]] constexpr auto segments() -> segments_type & {
    return m_segments;
  }

  [[nodiscard]] constexpr auto segments() const -> segments_type const & {
    return m_segments;
  }

  [[nodiscard]] auto split() const -> std::array<std::array<std::optional<exact_region>, 2>, 2> {
    if (m_opt_point1.has_value() && m_opt_point2.has_value()) {
      if (!is_point_close(*m_opt_point1, *m_opt_point2)) {
        auto r1 = m_split(*m_opt_point1);
        auto r2 = m_split(*m_opt_point2);
        if (r1[0].has_value() && r1[1].has_value() && r2[0].has_value() && r2[1].has_value()) {
          if (m_symmetric(*r1[0], *r2[1]) && m_symmetric(*r1[1], *r2[0])) {
            return {r1, {std::nullopt, std::nullopt}};
          }
        } else if (r1[0].has_value() && r2[1].has_value()) {
          if (m_symmetric(*r1[0], *r2[1])) {
            return {r1, {std::nullopt, std::nullopt}};
          }
        } else if (r1[1].has_value() && r2[0].has_value()) {
          if (m_symmetric(*r1[1], *r2[0])) {
            return {r1, {std::nullopt, std::nullopt}};
          }
        }
        return {r1, r2};
      } else {
        return {m_split(*m_opt_point1), {std::nullopt, std::nullopt}};
      }
    } else if (m_opt_point1.has_value()) {
      return {m_split(*m_opt_point1), {std::nullopt, std::nullopt}};
    } else {
      return {std::array<std::optional<exact_region>, 2>{std::nullopt, std::nullopt}, std::array<std::optional<exact_region>, 2>{std::nullopt, std::nullopt}};
    }
  }

 private:
  [[nodiscard]] auto m_split(point const &p) const -> std::array<std::optional<exact_region>, 2> {
    std::vector<segment> s1, s2;
    s1.reserve(m_segments.size());
    s2.reserve(m_segments.size());
    for (auto const &s : m_segments) {
      if (s.b.x < p.x || is_close(s.b.x, p.x)) {
        s1.push_back(s);
      } else if (s.a.x > p.x || is_close(s.a.x, p.x)) {
        s2.push_back(s);
      } else {
        s1.emplace_back(s.a, p);
        s2.emplace_back(p, s.b);
      }
    }
    auto ref1 = point{m_ref.x, p.y};
    auto ref2 = point{p.x, m_ref.y};

    if (!s1.empty() && !s2.empty()) {
      return {exact_region(std::move(s1), ref1), exact_region(std::move(s2), ref2)};
    } else if (!s1.empty()) {
      return {exact_region(std::move(s1), ref1), std::nullopt};
    } else if (!s2.empty()) {
      return {std::nullopt, exact_region(std::move(s2), ref2)};
    } else {
      return {std::nullopt, std::nullopt};
    }
  }

  [[nodiscard]] auto m_symmetric(exact_region const &r1, exact_region const &r2) const -> bool {
    if (r1.m_segments.size() != r2.m_segments.size()) {
      return false;
    }
    if (r1.empty()) {
      return true;
    }

    // Check if all segments have the same slope.
    auto iter1 = r1.m_segments.begin();
    auto last1 = r1.m_segments.end();
    auto iter2 = r2.m_segments.rbegin();
    for (; iter1 != last1; ++iter1, ++iter2) {
      auto ma = (iter1->b.y - iter1->a.y) / (iter1->b.x - iter1->a.x);
      auto mb = (iter2->a.x - iter2->b.x) / (iter2->a.y - iter2->b.y);
      if (!is_close(ma, mb)) {
        return false;
      }
    }

    // Finally check the segments distances
    iter1 = r1.m_segments.begin();
    iter2 = r2.m_segments.rbegin();
    auto p1 = point{r1.m_ref.x, iter1->a.y};
    auto p2 = point{iter2->b.x, r2.m_ref.y};
    for (; iter1 != last1; ++iter1, ++iter2) {
      // Same distance to previous point
      auto d1x = iter1->a.x - p1.x;
      auto d1y = p1.y - iter1->a.y;
      auto d2x = p2.x - iter2->b.x;
      auto d2y = iter2->b.y - p2.y;
      if (!is_close(d1x, d2y) || !is_close(d2x, d1y)) {
        return false;
      }

      // Same segment distance
      d1x = iter1->b.x - iter1->a.x;
      d1y = iter1->a.y - iter1->b.y;
      d2x = iter2->b.x - iter2->a.x;
      d2y = iter2->a.y - iter2->b.y;
      if (!is_close(d1x, d2y) || !is_close(d2x, d1y)) {
        return false;
      }

      p1 = iter1->b;
      p2 = iter2->a;
    }

    return true;
  }
};

class exact_model : public std::ranges::view_interface<exact_model> {
 private:
  std::vector<std::priority_queue<exact_region>> m_pqs;

 public:
  template <typename Segments>
  constexpr exact_model(Segments &&segments, point r)
      : m_pqs(1) {
    m_pqs.back().emplace(std::forward<Segments>(segments), r);
  }

  class sentinel {};

  class iterator {
   private:
    exact_model *parent = nullptr;
    double value;

   public:
    using iterator_category = std::input_iterator_tag;
    using value_type = double;
    using difference_type = std::ptrdiff_t;

    constexpr iterator() {}

    iterator(exact_model &p)
        : parent(&p)
        , value(parent->m_pqs.front().top().max_hv()) {}

    constexpr auto operator==(sentinel const &) const -> bool {
      return parent == NULL || parent->m_pqs.empty();
    }

    constexpr auto operator*() const -> value_type {
      return value;
    }

    auto operator++() -> iterator & {
      decltype(parent->m_pqs) pqs2;
      pqs2.reserve(parent->m_pqs.size());

      for (auto &pq : parent->m_pqs) {
        if (pq.top().max_hv() < value && !is_close(pq.top().max_hv(), value)) {
          continue;
        }

        // Takes only the lexicographically extreme points
        auto rs = pq.top().split();
        pq.pop();
        if ((rs[0][0].has_value() || rs[0][1].has_value()) && (rs[1][0].has_value() || rs[1][1].has_value())) {
          pqs2.push_back(pq);
          if (rs[0][0].has_value())
            pqs2.back().push(std::move(*rs[0][0]));
          if (rs[0][1].has_value())
            pqs2.back().push(std::move(*rs[0][1]));
          pqs2.push_back(std::move(pq));
          if (rs[1][0].has_value())
            pqs2.back().push(std::move(*rs[1][0]));
          if (rs[1][1].has_value())
            pqs2.back().push(std::move(*rs[1][1]));
        } else if ((rs[0][0].has_value() || rs[0][1].has_value())) {
          pqs2.push_back(std::move(pq));
          if (rs[0][0].has_value())
            pqs2.back().push(std::move(*rs[0][0]));
          if (rs[0][1].has_value())
            pqs2.back().push(std::move(*rs[0][1]));
        } else if ((rs[1][0].has_value() || rs[1][1].has_value())) {
          pqs2.push_back(std::move(pq));
          if (rs[1][0].has_value())
            pqs2.back().push(std::move(*rs[1][0]));
          if (rs[1][1].has_value())
            pqs2.back().push(std::move(*rs[1][1]));
        } else {
          if (!pq.empty())
            pqs2.push_back(std::move(pq));
        }
      }
      parent->m_pqs = std::move(pqs2);
      // std::cerr << parent->m_pqs.size() << "\n";
      value = std::ranges::max(std::ranges::views::transform(parent->m_pqs, [](auto const &q) { return q.top().max_hv(); }));
      return *this;
    }

    auto operator++(int) -> iterator & {
      return ++(*this);
    }
  };

  [[nodiscard]] auto begin() {
    return iterator(*this);
  }

  [[nodiscard]] constexpr auto end() {
    return sentinel{};
  }
};

class greedy_region {
 private:
  std::vector<segment> m_segments;
  point m_ref;
  std::optional<point> m_opt_point;
  double m_opt_point_hv;

 public:
  template <typename Segments>
  constexpr greedy_region(Segments &&segments, point ref)
      : m_segments(std::forward<Segments>(segments))
      , m_ref(ref)
      , m_opt_point(std::nullopt)
      , m_opt_point_hv(0.0) {
    for (auto const &s : m_segments) {
      auto const p = s.max_hv_point(m_ref);
      auto const hv = p.hv(m_ref);
      if (hv > m_opt_point_hv && !is_close(hv, m_opt_point_hv)) {
        m_opt_point_hv = hv;
        m_opt_point = p;
      }
    }
  }

  [[nodiscard]] constexpr auto operator<(greedy_region const &o) const -> bool {
    return m_opt_point_hv < o.m_opt_point_hv;
  }

  [[nodiscard]] constexpr auto max_hv() const -> double {
    return m_opt_point_hv;
  }

  [[nodiscard]] auto empty() const -> bool {
    return m_segments.empty();
  }

  [[nodiscard]] constexpr auto segments() -> std::vector<segment> & {
    return m_segments;
  }

  [[nodiscard]] constexpr auto segments() const -> std::vector<segment> const & {
    return m_segments;
  }

  [[nodiscard]] constexpr auto opt_point() const -> point const & {
    return *m_opt_point;
  }

  [[nodiscard]] auto split() const -> std::array<std::optional<greedy_region>, 2> {
    if (!m_opt_point.has_value())
      return {std::nullopt, std::nullopt};

    std::vector<segment> s1, s2;
    s1.reserve(m_segments.size());
    s2.reserve(m_segments.size());
    for (auto const &s : m_segments) {
      if (s.b.x < m_opt_point->x || is_close(s.b.x, m_opt_point->x)) {
        s1.push_back(s);
      } else if (s.a.x > m_opt_point->x || is_close(s.a.x, m_opt_point->x)) {
        s2.push_back(s);
      } else {
        s1.push_back({s.a, *m_opt_point});
        s2.push_back({*m_opt_point, s.b});
      }
    }
    auto ref1 = point{m_ref.x, m_opt_point->y};
    auto ref2 = point{m_opt_point->x, m_ref.y};
    if (!s1.empty() && !s2.empty()) {
      return {greedy_region(std::move(s1), ref1), greedy_region(std::move(s2), ref2)};
    } else if (!s1.empty()) {
      return {greedy_region(std::move(s1), ref1), std::nullopt};
    } else if (!s2.empty()) {
      return {std::nullopt, greedy_region(std::move(s2), ref2)};
    } else {
      return {std::nullopt, std::nullopt};
    }
  }
};

class greedy_model : public std::ranges::view_interface<greedy_model> {
 private:
  std::priority_queue<greedy_region> m_pq;

 public:
  template <typename Segments>
  constexpr greedy_model(Segments &&segments, point r) {
    m_pq.emplace(std::forward<Segments>(segments), r);
  }

  class sentinel {};

  class iterator {
   private:
    greedy_model *parent = nullptr;

   public:
    using iterator_category = std::input_iterator_tag;
    using value_type = std::pair<point, double>;
    using difference_type = std::ptrdiff_t;

    constexpr iterator() {}

    constexpr iterator(greedy_model &p)
        : parent(&p) {}

    constexpr auto operator==(sentinel const &) const -> bool {
      return parent == NULL || parent->m_pq.empty();
    }

    auto operator*() const -> value_type {
      return std::make_pair(parent->m_pq.top().opt_point(), parent->m_pq.top().max_hv());
    }

    auto operator++() -> iterator & {
      auto rs = parent->m_pq.top().split();
      parent->m_pq.pop();
      if (rs[0].has_value())
        parent->m_pq.push(std::move(*rs[0]));
      if (rs[1].has_value())
        parent->m_pq.push(std::move(*rs[1]));
      return *this;
    }

    auto operator++(int) -> iterator & {
      return ++(*this);
    }
  };

  [[nodiscard]] constexpr auto begin() {
    return iterator(*this);
  }

  [[nodiscard]] constexpr auto end() {
    return sentinel{};
  }
};

}  // namespace apm::details

namespace apm {
using details::exact_model;
using details::greedy_model;
using details::piecewise_segments;
using details::point;
using details::segment;
}  // namespace apm

#endif
