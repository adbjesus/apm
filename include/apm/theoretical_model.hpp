#ifndef APM_THEORETICAL_MODEL_HPP_
#define APM_THEORETICAL_MODEL_HPP_

#include "utils.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <concepts>
#include <iostream>
#include <numbers>
#include <optional>
#include <queue>
#include <ranges>
#include <type_traits>
#include <vector>

namespace apm::theoretical_model {

template <typename Scalar>
using point = std::array<Scalar, 2>;

template <typename Scalar>
using segment = std::array<point<Scalar>, 2>;

template <typename Scalar, Scalar atol = 1e-30, Scalar rtol = 1e-9>
struct is_close_fn {
  [[nodiscard]] constexpr auto operator()(Scalar const &lhs, Scalar const &rhs) const -> bool {
    return std::abs(lhs - rhs) <= atol + rtol * std::abs(rhs);
  }
};

template <typename Scalar, Scalar atol = 1e-30, Scalar rtol = 1e-9>
[[nodiscard]] constexpr auto is_close(Scalar const &lhs, Scalar const &rhs) -> bool {
  return std::abs(lhs - rhs) <= atol + rtol * std::abs(rhs);
}

/**
 * Returns the point that gives the maximal hypervolume contribution in
 * a segment `s` according to a reference point `r`.
 *
 * Assumptions:
 *   - s[0][0] < s[1][0]
 *   - s[0][1] > s[1][1]
 *   - s[0][0] >= r[0]
 *   - s[1][1] >= r[1]
 */
template <typename Scalar>
[[nodiscard]] constexpr auto segment_max_hv_point(segment<Scalar> const &s, point<Scalar> const &r) {
  assert(s[0][0] < s[1][0]);
  assert(s[0][1] > s[1][1]);
  assert(s[0][0] >= r[0]);
  assert(s[1][1] >= r[1]);

  auto const m = (s[1][1] - s[0][1]) / (s[1][0] - s[0][0]);
  auto const b = s[0][1] - m * s[0][0];
  auto const c0 = (r[1] - b) / m;
  auto const z0 = r[0] + (c0 - r[0]) / 2;
  if (z0 > s[1][0] || is_close(z0, s[1][0])) {
    return s[1];
  } else if (z0 < s[0][0] || is_close(z0, s[0][0])) {
    return s[0];
  } else {
    auto const c1 = m * r[0] + b;
    auto const z1 = r[1] + (c1 - r[1]) / 2;
    return point<Scalar>{z0, z1};
  }
}

/**
 * Returns the hypervolume of `p` with respect to `r`.
 *
 * Assumptions:
 *   - p[0] >= r[0]
 *   - p[1] >= r[1]
 */
template <typename Scalar>
[[nodiscard]] constexpr auto point_hv(point<Scalar> const &p, point<Scalar> const &r) {
  return (p[0] - r[0]) * (p[1] - r[1]);
}

/**
 * Generate `n` piecewise approximation segments according to parameter `d` and place them in
 * the range beggining at `iter`.
 */
template <typename Size, typename Scalar>
[[nodiscard]] constexpr auto piecewise_segments(Size n, Scalar d) -> std::vector<segment<Scalar>> {
  auto res = std::vector<segment<Scalar>>{};
  res.reserve(n);

  auto const pi2 = std::numbers::pi_v<Scalar> / static_cast<Scalar>(2);
  auto const p = static_cast<Scalar>(2) / d;
  auto const s = pi2 / static_cast<Scalar>(n);
  auto s0 = static_cast<Scalar>(0);
  auto s1 = static_cast<Scalar>(1);
  for (Size i = 1; i < n; ++i) {
    auto const te = std::min(s * static_cast<Scalar>(i), pi2);
    auto const e0 = std::pow(std::sin(te), p);
    auto const e1 = std::pow(std::cos(te), p);
    res.push_back({point<Scalar>{s0, s1}, point<Scalar>{e0, e1}});
    s0 = std::move(e0);
    s1 = std::move(e1);
  }
  res.push_back({point<Scalar>{s0, s1}, point<Scalar>{static_cast<Scalar>(1), static_cast<Scalar>(0)}});
  return res;
}

/**
 * Generate `n` segments according to parameter `d` and place them in
 * the range beggining at `iter`.
 */
template <typename Iter, typename Size, typename Scalar>
constexpr auto generate_segments(Iter iter, Size n, Scalar d) {
  auto const pi2 = std::numbers::pi_v<Scalar> / static_cast<Scalar>(2);
  auto const p = static_cast<Scalar>(2) / d;
  auto const s = pi2 / static_cast<Scalar>(n);
  auto s0 = static_cast<Scalar>(0);
  auto s1 = static_cast<Scalar>(1);
  for (Size i = 1; i < n; ++i, ++iter) {
    auto const te = std::min(s * static_cast<Scalar>(i), pi2);
    auto const e0 = std::pow(std::sin(te), p);
    auto const e1 = std::pow(std::cos(te), p);
    *iter = segment<Scalar>{point<Scalar>{s0, s1}, point<Scalar>{e0, e1}};
    s0 = std::move(e0);
    s1 = std::move(e1);
  }
  *iter = segment<Scalar>{point<Scalar>{s0, s1}, point<Scalar>{static_cast<Scalar>(1), static_cast<Scalar>(0)}};
}

template <typename Scalar = double, typename Segments = std::vector<segment<Scalar>>>
class region {
 public:
  using scalar_type = Scalar;
  using segments_type = Segments;
  using segment_type = typename segments_type::value_type;
  using point_type = typename segment_type::value_type;

  template <typename SegmentsRange, typename P>
  constexpr region(SegmentsRange &&segments, P &&ref)
      : m_segments(std::forward<SegmentsRange>(segments))
      , m_ref(std::forward<P>(ref))
      , m_p1(std::nullopt)
      , m_p2(std::nullopt)
      , m_hv(0) {
    if (m_segments.empty())
      return;
    m_p1 = segment_max_hv_point(m_segments.front(), m_ref);
    m_hv = point_hv(*m_p1, m_ref);
    for (auto it = std::next(m_segments.begin()); it != m_segments.end(); ++it) {
      auto const p = segment_max_hv_point(*it, m_ref);
      auto const hv = point_hv(p, m_ref);
      if (is_close(hv, m_hv)) {
        m_p2 = std::move(p);
      } else if (hv > m_hv) {
        m_hv = hv;
        m_p1 = std::move(p);
        m_p2 = std::nullopt;
      }
    }
  }

  [[nodiscard]] constexpr auto operator<(region const &o) const -> bool {
    return m_hv < o.m_hv;
  }

  [[nodiscard]] constexpr auto hv() const {
    return m_hv;
  }

  [[nodiscard]] constexpr auto empty() const -> bool {
    return m_segments.emtpy();
  }

  [[nodiscard]] constexpr auto segments() -> segments_type & {
    return m_segments;
  }

  [[nodiscard]] constexpr auto segments() const -> segments_type const & {
    return m_segments;
  }

  [[nodiscard]] constexpr auto split() const -> std::array<std::array<std::optional<region>, 2>, 2> {
    if (m_p1.has_value() && m_p2.has_value()) {
      if (!is_close((*m_p1)[0], (*m_p2)[0]) || !is_close((*m_p1)[1], (*m_p2)[1])) {
        return {m_split(*m_p1), m_split(*m_p2)};
      } else {
        return {m_split(*m_p1), {std::nullopt, std::nullopt}};
      }
    } else if (m_p1.has_value()) {
      return {m_split(*m_p1), {std::nullopt, std::nullopt}};
    } else {
      return {std::array<std::optional<region>, 2>{std::nullopt, std::nullopt}, std::array<std::optional<region>, 2>{std::nullopt, std::nullopt}};
    }
  }

 private:
  [[nodiscard]] constexpr auto m_split(point_type const &p) const -> std::array<std::optional<region>, 2> {
    segments_type s1, s2;
    s1.reserve(m_segments.size());
    s2.reserve(m_segments.size());
    for (auto const &s : m_segments) {
      if (s[1][0] < p[0] || is_close(s[1][0], p[0])) {
        s1.push_back(s);
      } else if (s[0][0] > p[0] || is_close(s[0][0], p[0])) {
        s2.push_back(s);
      } else {
        s1.push_back({s[0], p});
        s2.push_back({p, s[1]});
      }
    }
    auto ref1 = point_type{m_ref[0], p[1]};
    auto ref2 = point_type{p[0], m_ref[1]};

    if (!s1.empty() && !s2.empty()) {
      return {region(std::move(s1), ref1), region(std::move(s2), ref2)};
    } else if (!s1.empty()) {
      return {region(std::move(s1), ref1), std::nullopt};
    } else if (!s2.empty()) {
      return {std::nullopt, region(std::move(s2), ref2)};
    } else {
      return {std::nullopt, std::nullopt};
    }
  }

  segments_type m_segments;
  point_type m_ref;
  std::optional<point_type> m_p1;
  std::optional<point_type> m_p2;
  scalar_type m_hv;
};

template <typename Scalar>
class exact_model : public std::ranges::view_interface<exact_model<Scalar>> {
 public:
  using scalar_type = Scalar;
  using point_type = point<scalar_type>;

 private:
  using segments_type = std::vector<segment<scalar_type>>;
  using region_type = region<scalar_type, segments_type>;
  using priority_queue_type = std::priority_queue<region_type>;

  std::vector<priority_queue_type> m_pqs;

 public:
  constexpr exact_model(segments_type &&segments, point_type r)
      : m_pqs(1) {
    m_pqs.back().emplace(std::move(segments), r);
  }

  constexpr exact_model(segments_type const &segments, point_type r)
      : m_pqs(1) {
    m_pqs.back().emplace(segments, r);
  }

  template <typename... Args>
  constexpr exact_model(Args &&...args, point_type r)
      : m_pqs(1) {
    m_pqs.back().emplace(segments_type(std::forward<Args>(args)...), r);
  }

  class sentinel {};

  class iterator {
   private:
    exact_model *parent = nullptr;
    scalar_type value;

   public:
    using iterator_category = std::input_iterator_tag;
    using value_type = scalar_type;
    using difference_type = std::ptrdiff_t;

    constexpr iterator() {}

    constexpr iterator(exact_model &p)
        : parent(&p)
        , value(parent->m_pqs.front().top().hv()) {}

    constexpr auto operator==(sentinel const &) const -> bool {
      return parent == NULL || parent->m_pqs.empty();
    }

    constexpr auto operator*() const -> value_type {
      return value;
    }

    constexpr auto operator++() -> iterator & {
      std::vector<priority_queue_type> pqs2;
      pqs2.reserve(parent->m_pqs.size());

      for (auto &pq : parent->m_pqs) {
        if (pq.top().hv() < value || !is_close(pq.top().hv(), value)) {
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
      value = std::ranges::max(std::ranges::views::transform(parent->m_pqs, [](auto const &q) { return q.top().hv(); }));
      return *this;
    }

    constexpr auto operator++(int) -> iterator & {
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

// template <typename Iter, typename Segments, typename Size, typename Scalar>
// constexpr auto algorithmic_model(Iter iter, Size n, Segments &&segments, Scalar r0, Scalar r1) {
//   using scalar_type = typename std::remove_cvref_t<Scalar>;
//   using segments_type = typename std::remove_cvref_t<Segments>;
//   using segment_type = typename segments_type::value_type;
//   using point_type = typename segment_type::value_type;
//   using region_type = region<scalar_type, segments_type>;

//   std::priority_queue<region_type> regions;
//   regions.emplace(std::forward<Segments>(segments), point_type{r0, r1});

//   auto pqs = std::vector{std::move(regions)};
//   for (Size i = 0; i < n; ++i, ++iter) {
//     auto hv = std::ranges::max(std::ranges::views::transform(pqs, [](auto const &q) { return q.top().hv(); }));
//     *iter = hv;

//     decltype(pqs) pqs2;
//     pqs2.reserve(pqs.size());

//     for (auto &pq : pqs) {
//       if (pq.top().hv() < hv || !is_close(pq.top().hv(), hv)) {
//         continue;
//       }

//       // Takes only the lexicographically extreme points
//       auto rs = pq.top().split();
//       pq.pop();
//       if ((rs[0][0].has_value() || rs[0][1].has_value()) && (rs[1][0].has_value() || rs[1][1].has_value())) {
//         pqs2.push_back(pq);
//         if (rs[0][0].has_value())
//           pqs2.back().push(std::move(*rs[0][0]));
//         if (rs[0][1].has_value())
//           pqs2.back().push(std::move(*rs[0][1]));
//         pqs2.push_back(std::move(pq));
//         if (rs[1][0].has_value())
//           pqs2.back().push(std::move(*rs[1][0]));
//         if (rs[1][1].has_value())
//           pqs2.back().push(std::move(*rs[1][1]));
//       } else if ((rs[0][0].has_value() || rs[0][1].has_value())) {
//         pqs2.push_back(std::move(pq));
//         if (rs[0][0].has_value())
//           pqs2.back().push(std::move(*rs[0][0]));
//         if (rs[0][1].has_value())
//           pqs2.back().push(std::move(*rs[0][1]));
//       } else if ((rs[1][0].has_value() || rs[1][1].has_value())) {
//         pqs2.push_back(std::move(pq));
//         if (rs[1][0].has_value())
//           pqs2.back().push(std::move(*rs[1][0]));
//         if (rs[1][1].has_value())
//           pqs2.back().push(std::move(*rs[1][1]));
//       } else {
//         if (!pq.empty())
//           pqs2.push_back(std::move(pq));
//       }
//     }
//     pqs = std::move(pqs2);
//   }
// }

template <typename Scalar, typename Segments>
class region_greedy {
 public:
  using scalar_type = Scalar;
  using segments_type = Segments;
  using segment_type = typename segments_type::value_type;
  using point_type = typename segment_type::value_type;

  template <typename S, typename P>
  constexpr region_greedy(S &&segments, P &&ref)
      : m_segments(std::forward<S>(segments))
      , m_ref(std::forward<P>(ref))
      , m_p(std::nullopt)
      , m_hv(static_cast<Scalar>(0)) {
    for (auto const &s : m_segments) {
      auto const p = segment_max_hv_point(s, m_ref);
      auto const hv = (p[0] - m_ref[0]) * (p[1] - m_ref[1]);
      if (hv > m_hv && !is_close(hv, m_hv)) {
        m_hv = hv;
        m_p = p;
      }
    }
  }

  [[nodiscard]] constexpr auto operator<(region_greedy const &o) const -> bool {
    return m_hv < o.m_hv;
  }

  [[nodiscard]] constexpr auto hv() const -> scalar_type const & {
    return m_hv;
  }

  [[nodiscard]] constexpr auto empty() const {
    return m_segments.emtpy();
  }

  [[nodiscard]] constexpr auto segments() -> segments_type & {
    return m_segments;
  }

  [[nodiscard]] constexpr auto segments() const -> segments_type const & {
    return m_segments;
  }

  [[nodiscard]] constexpr auto opt_point() const -> point_type const & {
    return *m_p;
  }

  [[nodiscard]] constexpr auto split() const -> std::array<std::optional<region_greedy>, 2> {
    if (!m_p.has_value())
      return {std::nullopt, std::nullopt};

    segments_type s1, s2;
    s1.reserve(m_segments.size());
    s2.reserve(m_segments.size());
    for (auto const &s : m_segments) {
      if (s[1][0] < (*m_p)[0] || is_close(s[1][0], (*m_p)[0])) {
        s1.push_back(s);
      } else if (s[0][0] > (*m_p)[0] || is_close(s[0][0], (*m_p)[0])) {
        s2.push_back(s);
      } else {
        s1.push_back({s[0], *m_p});
        s2.push_back({*m_p, s[1]});
      }
    }
    auto ref1 = point_type{m_ref[0], (*m_p)[1]};
    auto ref2 = point_type{(*m_p)[0], m_ref[1]};
    if (!s1.empty() && !s2.empty()) {
      return {region_greedy(std::move(s1), ref1), region_greedy(std::move(s2), ref2)};
    } else if (!s1.empty()) {
      return {region_greedy(std::move(s1), ref1), std::nullopt};
    } else if (!s2.empty()) {
      return {std::nullopt, region_greedy(std::move(s2), ref2)};
    } else {
      return {std::nullopt, std::nullopt};
    }
  }

 private:
  segments_type m_segments;
  point_type m_ref;
  std::optional<point_type> m_p;
  scalar_type m_hv;
};

template <typename Scalar>
class greedy_model : public std::ranges::view_interface<greedy_model<Scalar>> {
 public:
  using scalar_type = Scalar;
  using point_type = point<scalar_type>;

 private:
  using segments_type = std::vector<segment<scalar_type>>;
  using region_type = region_greedy<scalar_type, segments_type>;
  using priority_queue_type = std::priority_queue<region_type>;

  priority_queue_type m_pq;

 public:
  constexpr greedy_model(segments_type &&segments, point_type r) {
    m_pq.emplace(std::move(segments), r);
  }

  constexpr greedy_model(segments_type const &segments, point_type r) {
    m_pq.emplace(segments, r);
  }

  template <typename... Args>
  constexpr greedy_model(Args &&...args, point_type r) {
    m_pq.emplace(segments_type(std::forward<Args>(args)...), r);
  }

  class sentinel {};

  class iterator {
   private:
    greedy_model *parent = nullptr;

   public:
    using iterator_category = std::input_iterator_tag;
    using value_type = std::pair<point<scalar_type>, scalar_type>;
    using difference_type = std::ptrdiff_t;

    constexpr iterator() {}

    constexpr iterator(greedy_model &p)
        : parent(&p) {}

    constexpr auto operator==(sentinel const &) const -> bool {
      return parent == NULL || parent->m_pq.empty();
    }

    constexpr auto operator*() const -> value_type {
      return make_pair(parent->m_pq.top().opt_point(), parent->m_pq.top().hv());
    }

    constexpr auto operator++() -> iterator & {
      auto rs = parent->m_pq.top().split();
      parent->m_pq.pop();
      if (rs[0].has_value())
        parent->m_pq.push(std::move(*rs[0]));
      if (rs[1].has_value())
        parent->m_pq.push(std::move(*rs[1]));
      return *this;
    }

    constexpr auto operator++(int) -> iterator & {
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

// template <std::floating_point T>
// std::vector<point<T>> tapm_greedy_points(size_t n, size_t l, T d, T rx, T ry, int ulp = 2) {
//   auto result = std::vector<point<T>>();
//   result.reserve(n);

//   auto segments = generate_segments(l, d);

//   std::priority_queue<region_greedy<T>> pq;
//   pq.emplace(std::move(segments), point<T>(rx, ry), ulp);
//   for (size_t i = 0; i < n && !pq.empty(); ++i) {
//     auto const &r = pq.top();
//     result.push_back(r.opt());
//     pq.pop();
//     auto [r0, r1] = r.split();
//     if (r0.has_value())
//       pq.push(std::move(*r0));
//     if (r1.has_value())
//       pq.push(std::move(*r1));
//   }

//   return result;
// }

}  // namespace apm::theoretical_model

namespace apm {
using theoretical_model::exact_model;
using theoretical_model::greedy_model;
using theoretical_model::piecewise_segments;
using theoretical_model::point;
using theoretical_model::segment;

// using theoretical_model::analytical_model;
}  // namespace apm

#endif
