#ifndef APM_APM_HPP_
#define APM_APM_HPP_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <concepts>
#include <numbers>
#include <optional>
#include <queue>
#include <type_traits>
#include <vector>

#include "utils.hpp"

namespace apm {

template <std::floating_point T>
struct point {
  T x;
  T y;

  point(T x, T y)
      : x(x)
      , y(y) {}
};

// Assumptions: start.x < end.x
template <std::floating_point T>
struct segment {
 public:
  using float_type = T;
  using point_type = point<T>;

  segment(point_type s, point_type e)
      : s(s)
      , e(e)
      , m((e.y - s.y) / (e.x - s.x))
      , b(s.y - m * s.x) {}

  [[nodiscard]] constexpr auto max_hv_point(point_type const &r, int ulp) const -> point_type {
    float_type c0 = (r.y - b) / m;
    float_type z0 = r.x + (c0 - r.x) / 2;
    if (flt_ge(z0, e.x, ulp)) {
      return e;
    } else if (flt_le(z0, s.x, ulp)) {
      return s;
    } else {
      float_type c1 = m * r.x + b;
      float_type z1 = r.y + (c1 - r.y) / 2;
      return point_type(z0, z1);
    }
  }

 public:
  point_type s;
  point_type e;

 private:
  float_type m;
  float_type b;
};

template <std::floating_point T>
[[nodiscard]] constexpr auto generate_segments(size_t n, T d) -> std::vector<segment<T>> {
  std::vector<segment<T>> segments;
  segments.reserve(n);
  const auto pi2 = std::numbers::pi_v<T> / T(2);
  const auto p = T(2) / d;
  const auto s = pi2 / T(n);
  auto s0 = T(0);
  auto s1 = T(1);
  for (size_t i = 1; i < n; ++i) {
    const auto te = std::min(s * T(i), pi2);
    const auto e0 = std::pow(std::sin(te), p);
    const auto e1 = std::pow(std::cos(te), p);
    segments.push_back(segment<T>{{s0, s1}, {e0, e1}});
    s0 = e0;
    s1 = e1;
  }
  segments.push_back(segment<T>{{s0, s1}, {1, 0}});
  return segments;
}

template <std::floating_point T>
class region {
 public:
  using float_type = T;
  using point_type = point<float_type>;
  using segment_type = segment<float_type>;
  using segments_type = std::vector<segment_type>;

  constexpr region(segments_type &&segments, point_type ref, int ulp)
      : m_segments(std::move(segments))
      , m_opt1(std::nullopt)
      , m_opt2(std::nullopt)
      , m_ref(ref)
      , m_hv(0)
      , m_ulp(ulp) {
    if (m_segments.empty())
      return;
    m_opt1 = m_segments.begin()->max_hv_point(m_ref, m_ulp);
    m_hv = (m_opt1->x - m_ref.x) * (m_opt1->y - m_ref.y);
    for (auto it = std::next(m_segments.begin()); it != m_segments.end(); ++it) {
      auto z = it->max_hv_point(m_ref, m_ulp);
      float_type hv = (z.x - m_ref.x) * (z.y - m_ref.y);
      if (flt_eq(hv, m_hv, m_ulp)) {
        if (!flt_eq(m_opt1->x, z.x, ulp) || !flt_eq(m_opt1->y, z.y, ulp)) {
          m_opt2 = std::move(z);
        }
      } else if (hv > m_hv) {
        m_hv = hv;
        m_opt1 = std::move(z);
        m_opt2 = std::nullopt;
      }
    }
  }

  [[nodiscard]] constexpr auto operator<(region<T> const &o) const -> bool {
    return flt_lt(m_hv, o.m_hv, m_ulp);
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

  [[nodiscard]] constexpr auto split() const -> std::array<std::array<std::optional<region<float_type>>, 2>, 2> {
    return {m_split(m_opt1), m_split(m_opt2)};
  }

 private:
  [[nodiscard]] constexpr auto m_split(std::optional<point_type> const &p) const
      -> std::array<std::optional<region<float_type>>, 2> {
    if (!p.has_value())
      return {std::nullopt, std::nullopt};

    std::vector<segment_type> s1;
    std::vector<segment_type> s2;
    for (auto const &s : m_segments) {
      if (flt_le(s.e.x, p->x, m_ulp)) {
        s1.push_back(s);
      } else if (flt_ge(s.s.x, p->x, m_ulp)) {
        s2.push_back(s);
      } else {
        s1.push_back(segment{s.s, *p});
        s2.push_back(segment{*p, s.e});
      }
    }
    auto ref1 = point_type{m_ref.x, p->y};
    auto ref2 = point_type{p->x, m_ref.y};

    if (!s1.empty() && !s2.empty()) {
      return {region<float_type>(std::move(s1), ref1, m_ulp), region<float_type>(std::move(s2), ref2, m_ulp)};
    } else if (!s1.empty()) {
      return {region<float_type>(std::move(s1), ref1, m_ulp), std::nullopt};
    } else if (!s2.empty()) {
      return {std::nullopt, region<float_type>(std::move(s2), ref2, m_ulp)};
    } else {
      return {std::nullopt, std::nullopt};
    }
  }

  std::vector<segment_type> m_segments;
  std::optional<point_type> m_opt1;
  std::optional<point_type> m_opt2;
  point_type m_ref;
  float_type m_hv;
  int m_ulp;
};

template <std::floating_point T>
std::vector<T> tapm(size_t n, size_t l, T d, T rx, T ry, int ulp = 2) {
  auto result = std::vector<T>();
  result.reserve(n);

  auto segments = generate_segments(l, d);

  std::priority_queue<region<T>> regions;
  regions.emplace(std::move(segments), point<T>(rx, ry), ulp);

  std::vector<std::priority_queue<region<T>>> pqs;
  pqs.push_back(std::move(regions));
  for (size_t i = 0; i < n; ++i) {
    T hv = 0;
    for (auto const &pq : pqs) {
      hv = std::max(hv, pq.top().hv());
    }
    result.push_back(hv);
    std::vector<std::priority_queue<region<T>>> pqs2;
    pqs2.reserve(pqs.size());
    for (auto &pq : pqs) {
      if (flt_lt(pq.top().hv(), hv, ulp)) {
        continue;
      }
      auto r = std::move(const_cast<region<T> &>(pq.top()));
      pq.pop();
      // Takes only the lexicographically extreme points
      auto rs = r.split();
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
    pqs = std::move(pqs2);
  }

  return result;
}

template <std::floating_point T>
class region_greedy {
 public:
  using float_type = T;
  using point_type = point<float_type>;
  using segment_type = segment<float_type>;
  using segments_type = std::vector<segment_type>;

  constexpr region_greedy(std::vector<segment_type> &&segments, point_type ref, int ulp)
      : m_segments(std::move(segments))
      , m_opt(std::nullopt)
      , m_ref(ref)
      , m_hv(0)
      , m_ulp(ulp) {
    for (auto const &s : m_segments) {
      auto z = s.max_hv_point(m_ref, m_ulp);
      float_type hv = (z.x - m_ref.x) * (z.y - m_ref.y);
      if (flt_gt(hv, m_hv, ulp)) {
        m_hv = hv;
        m_opt = z;
      }
    }
  }

  [[nodiscard]] constexpr auto operator<(region_greedy<float_type> const &o) const -> bool {
    return flt_lt(m_hv, o.m_hv, m_ulp);
  }

  [[nodiscard]] constexpr auto hv() const {
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

  [[nodiscard]] constexpr auto opt() const -> point_type {
    return *m_opt;
  }

  [[nodiscard]] constexpr auto split() const -> std::array<std::optional<region_greedy<float_type>>, 2> {
    if (!m_opt.has_value())
      return {std::nullopt, std::nullopt};

    std::vector<segment_type> s1;
    std::vector<segment_type> s2;
    for (auto const &s : m_segments) {
      if (flt_le(s.e.x, m_opt->x, m_ulp)) {
        s1.push_back(s);
      } else if (flt_ge(s.s.x, m_opt->x, m_ulp)) {
        s2.push_back(s);
      } else {
        s1.push_back(segment{s.s, *m_opt});
        s2.push_back(segment{*m_opt, s.e});
      }
    }
    auto ref1 = point_type(m_ref.x, m_opt->y);
    auto ref2 = point_type(m_opt->x, m_ref.y);
    if (!s1.empty() && !s2.empty()) {
      return {region_greedy<float_type>(std::move(s1), ref1, m_ulp),
              region_greedy<float_type>(std::move(s2), ref2, m_ulp)};
    } else if (!s1.empty()) {
      return {region_greedy<float_type>(std::move(s1), ref1, m_ulp), std::nullopt};
    } else if (!s2.empty()) {
      return {std::nullopt, region_greedy<float_type>(std::move(s2), ref2, m_ulp)};
    } else {
      return {std::nullopt, std::nullopt};
    }
  }

 private:
  segments_type m_segments;
  std::optional<point_type> m_opt;
  point_type m_ref;
  float_type m_hv;
  int m_ulp;
};

template <std::floating_point T>
std::vector<T> tapm_greedy(size_t n, size_t l, T d, T rx, T ry, int ulp = 2) {
  auto result = std::vector<T>();
  result.reserve(n);

  auto segments = generate_segments(l, d);

  std::priority_queue<region_greedy<T>> pq;
  pq.emplace(std::move(segments), point<T>(rx, ry), ulp);
  for (size_t i = 0; i < n && !pq.empty(); ++i) {
    auto r = pq.top();
    pq.pop();
    result.push_back(r.hv());
    auto [r0, r1] = r.split();
    if (r0.has_value())
      pq.push(std::move(*r0));
    if (r1.has_value())
      pq.push(std::move(*r1));
  }

  return result;
}

template <std::floating_point T>
std::vector<point<T>> tapm_greedy_points(size_t n, size_t l, T d, T rx, T ry, int ulp = 2) {
  auto result = std::vector<point<T>>();
  result.reserve(n);

  auto segments = generate_segments(l, d);

  std::priority_queue<region_greedy<T>> pq;
  pq.emplace(std::move(segments), point<T>(rx, ry), ulp);
  for (size_t i = 0; i < n && !pq.empty(); ++i) {
    auto const &r = pq.top();
    result.push_back(r.opt());
    pq.pop();
    auto [r0, r1] = r.split();
    if (r0.has_value())
      pq.push(std::move(*r0));
    if (r1.has_value())
      pq.push(std::move(*r1));
  }

  return result;
}

}  // namespace apm

#endif
