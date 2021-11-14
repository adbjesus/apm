#include <catch2/catch.hpp>

#include <apm/apm.hpp>

TEST_CASE("is_close to zero", "[apm::details][is_close]") {
  double rtol = 0.0;
  double atol = GENERATE(1e-40, 1e-35, 1e-30);
  REQUIRE(apm::details::is_close(atol, 0.0, rtol, atol) == true);
  REQUIRE(apm::details::is_close(atol * 2.0, 0.0, rtol, atol) == false);
}

TEST_CASE("is_close", "[apm::details][is_close]") {
  double atol = 0.0;
  double rtol = GENERATE(1e-8, 1e-6, 1e-4);
  double val = GENERATE(1e-20, 1e-15, 1e-10, 1e-5, 1e-1, 1);
  REQUIRE(apm::details::is_close(val, val, rtol, atol) == true);
  REQUIRE(apm::details::is_close(val, val + val * rtol, rtol, atol) == true);
  REQUIRE(apm::details::is_close(val + val * rtol, val, rtol, atol) == true);
  REQUIRE(apm::details::is_close(val, val + val * rtol * 1.1, rtol, atol) == false);
  REQUIRE(apm::details::is_close(val + val * rtol * 1.1, val, rtol, atol) == false);
}
