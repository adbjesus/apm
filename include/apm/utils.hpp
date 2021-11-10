#pragma once

#include <cmath>
#include <concepts>
#include <numeric>

namespace apm {

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

} // namespace apm
