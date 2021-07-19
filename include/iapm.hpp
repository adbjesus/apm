#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <concepts>
#include <iostream>
#include <numbers>
#include <queue>
#include <type_traits>
#include <vector>

namespace iapm {
// taken from https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
template <std::floating_point T> bool flt_eq(T x, T y, int ulp) {
  // the machine epsilon has to be scaled to the magnitude of the values used
  // and multiplied by the desired precision in ULPs (units in the last place)
  return std::abs(x - y) <
             std::numeric_limits<T>::epsilon() * std::abs(x + y) * ulp
         // unless the result is subnormal
         || std::abs(x - y) < std::numeric_limits<T>::min();
}

template <std::floating_point T> bool flt_lt(T x, T y, int ulp) {
  return !flt_eq(x, y, ulp) && x < y;
}

template <std::floating_point T> bool flt_le(T x, T y, int ulp) {
  return flt_eq(x, y, ulp) || x < y;
}

template <std::floating_point T> bool flt_gt(T x, T y, int ulp) {
  return !flt_eq(x, y, ulp) && x > y;
}

template <std::floating_point T> bool flt_ge(T x, T y, int ulp) {
  return flt_eq(x, y, ulp) || x > y;
}

template <std::floating_point T> struct segment {
  std::array<T, 2> s;
  std::array<T, 2> e;
};

template <std::floating_point T> class region {
public:
  using float_type = T;
  using point_type = std::array<float_type, 2>;
  using segment_type = segment<float_type>;

  region(std::vector<segment_type> &&segments, point_type ref, int ulp)
      : m_segments(std::move(segments)), m_ref((std::move(ref))), m_ulp(ulp) {
    m_hv = -1;
    for (auto const &s : m_segments) {
      // Compute optimal point for segment
      assert(s.e[0] > s.s[0] && s.e[1] < s.s[1]);
      float_type m = (s.e[1] - s.s[1]) / (s.e[0] - s.s[0]);
      float_type b = s.s[1] - m * s.s[0];
      float_type c0 = (m_ref[1] - b) / m;
      float_type c1 = m * m_ref[0] + b;
      float_type z0 = m_ref[0] + (c0 - m_ref[0]) / 2;
      float_type z1 = m_ref[1] + (c1 - m_ref[1]) / 2;
      point_type z;
      if (flt_ge(z0, s.e[0], m_ulp)) {
        z = s.e;
      } else if (z0 < s.s[0]) {
        z = s.s;
      } else {
        z = point_type{z0, z1};
      }
      float_type hv = (z[0] - m_ref[0]) * (z[1] - m_ref[1]);
      if (flt_eq(hv, m_hv, m_ulp)) {
        bool equal = false;
        for (auto const &p : m_opt) {
          if (flt_eq(p[0], z[0], ulp) && flt_eq(p[1], z[1], ulp)) {
            equal = true;
            break;
          }
        }
        if (!equal) {
          m_opt.push_back(z);
        }
      } else if (hv > m_hv) {
        m_hv = hv;
        m_opt.clear();
        m_opt.push_back(z);
      }
    }
  }

  bool operator<(region<T> const &o) const {
    return flt_lt(m_hv, o.m_hv, m_ulp);
  }

  auto hv() const { return m_hv; }

  auto begin() const { return m_opt.begin(); }

  auto end() const { return m_opt.end(); }

  auto size() const { return m_opt.size(); }

  std::vector<region<float_type>> split(point_type const &p) const {
    std::vector<segment_type> s1;
    std::vector<segment_type> s2;
    for (auto const &s : m_segments) {
      if (flt_le(s.e[0], p[0], m_ulp)) {
        s1.push_back(s);
      } else if (flt_ge(s.s[0], p[0], m_ulp)) {
        s2.push_back(s);
      } else {
        s1.push_back(segment{s.s, p});
        s2.push_back(segment{p, s.e});
      }
    }
    std::vector<region<float_type>> ret;
    if (!s1.empty())
      ret.emplace_back(std::move(s1), point_type{m_ref[0], p[1]}, m_ulp);
    if (!s2.empty())
      ret.emplace_back(std::move(s2), point_type{p[0], m_ref[1]}, m_ulp);
    return ret;
  }

private:
  std::vector<segment_type> m_segments;
  std::vector<point_type> m_opt;
  point_type m_ref;
  float_type m_hv;
  int m_ulp;
};

template <std::floating_point T>
std::vector<segment<T>> generate_segments(size_t n, T d) {
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
std::vector<T> iapm(size_t n, size_t l, T d, std::array<T, 2> r, int ulp = 2) {
  auto result = std::vector<T>();
  result.reserve(n);

  auto segments = generate_segments(l, d);

  std::priority_queue<region<T>> regions;
  regions.emplace(std::move(segments), r, ulp);

  std::vector<std::priority_queue<region<T>>> pqs;
  pqs.push_back(std::move(regions));
  for (size_t i = 0; i < n; ++i) {
    T hv = 0;
    for (auto &pq : pqs) {
      hv = std::max(hv, pq.top().hv());
    }
    result.push_back(hv);
    pqs.erase(std::remove_if(pqs.begin(), pqs.end(),
                             [hv, ulp](auto const &pq) {
                               return flt_lt(pq.top().hv(), hv, ulp);
                             }),
              pqs.end());
    std::cerr << i << " " << pqs.size() << std::endl;
    std::vector<std::priority_queue<region<T>>> pqs2;
    for (auto &pq : pqs) {
      auto r = pq.top();
      pq.pop();
      for (auto const &p : r) {
        pqs2.push_back(pq);
        for (auto aux : r.split(p)) {
          pqs2.back().push(std::move(aux));
        }
      }
    }
    pqs = std::move(pqs2);
  }

  return result;
}

template <std::floating_point T> class region_greedy {
public:
  using float_type = T;
  using point_type = std::array<float_type, 2>;
  using segment_type = segment<float_type>;

  region_greedy(std::vector<segment_type> &&segments, point_type ref, int ulp)
      : m_segments(std::move(segments)), m_ref((std::move(ref))), m_ulp(ulp) {
    m_hv = -1;
    for (auto const &s : m_segments) {
      // Compute optimal point for segment
      assert(s.e[0] > s.s[0] && s.e[1] < s.s[1]);
      float_type m = (s.e[1] - s.s[1]) / (s.e[0] - s.s[0]);
      float_type b = s.s[1] - m * s.s[0];
      float_type c0 = (m_ref[1] - b) / m;
      float_type c1 = m * m_ref[0] + b;
      float_type z0 = m_ref[0] + (c0 - m_ref[0]) / 2;
      float_type z1 = m_ref[1] + (c1 - m_ref[1]) / 2;
      point_type z;
      if (flt_ge(z0, s.e[0], m_ulp)) {
        z = s.e;
      } else if (z0 < s.s[0]) {
        z = s.s;
      } else {
        z = point_type{z0, z1};
      }
      float_type hv = (z[0] - m_ref[0]) * (z[1] - m_ref[1]);
      if (flt_gt(hv, m_hv, ulp)) {
        m_hv = hv;
        m_opt = z;
      }
    }
  }

  bool operator<(region_greedy<T> const &o) const {
    return flt_lt(m_hv, o.m_hv, m_ulp);
  }

  auto hv() const { return m_hv; }

  std::vector<region_greedy<float_type>> split() const {
    auto const &p = m_opt;
    std::vector<segment_type> s1;
    std::vector<segment_type> s2;
    for (auto const &s : m_segments) {
      if (flt_le(s.e[0], p[0], m_ulp)) {
        s1.push_back(s);
      } else if (flt_ge(s.s[0], p[0], m_ulp)) {
        s2.push_back(s);
      } else {
        s1.push_back(segment{s.s, p});
        s2.push_back(segment{p, s.e});
      }
    }
    std::vector<region_greedy<float_type>> ret;
    if (!s1.empty())
      ret.emplace_back(std::move(s1), point_type{m_ref[0], p[1]}, m_ulp);
    if (!s2.empty())
      ret.emplace_back(std::move(s2), point_type{p[0], m_ref[1]}, m_ulp);
    return ret;
  }

private:
  std::vector<segment_type> m_segments;
  point_type m_opt;
  point_type m_ref;
  float_type m_hv;
  int m_ulp;
};

template <std::floating_point T>
std::vector<T> iapm_greedy(size_t n, size_t l, T d, std::array<T, 2> r,
                           int ulp = 2) {
  auto result = std::vector<T>();
  result.reserve(n);

  auto segments = generate_segments(l, d);

  std::priority_queue<region_greedy<T>> pq;
  pq.emplace(std::move(segments), r, ulp);
  for (size_t i = 0; i < n; ++i) {
    auto r = pq.top();
    pq.pop();
    result.push_back(r.hv());
    for (auto aux : r.split()) {
      pq.push(std::move(aux));
    }
  }

  return result;
}
//
} // namespace iapm
