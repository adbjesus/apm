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

TEST_CASE("piecewise_segments", "[apm][piecewise_segments]") {
  size_t n = GENERATE(static_cast<size_t>(1), 2, 3, 4, 5, 10, 100, 1000);
  double d = GENERATE(0.1, 0.5, 1.0, 2.0, 10.0);
  auto segments = apm::piecewise_segments(n, d);
  // Correct size
  REQUIRE(segments.size() == n);
  if (n == 1) {
    REQUIRE(segments[0].a.x == 0.0);
    REQUIRE(segments[0].a.y == 1.0);
    REQUIRE(segments[0].b.x == 1.0);
    REQUIRE(segments[0].b.y == 0.0);
  } else if (n == 2) {
    auto p = std::pow(std::sin(std::numbers::pi_v<double> / 4.0), 2.0 / d);
    REQUIRE(segments[0].a.x == 0.0);
    REQUIRE(segments[0].a.y == 1.0);
    INFO(p);
    INFO(segments[0].b.x);
    REQUIRE(apm::details::is_close(segments[0].b.x, p));
    REQUIRE(apm::details::is_close(segments[0].b.y, p));
    REQUIRE(apm::details::is_close(segments[1].a.x, p));
    REQUIRE(apm::details::is_close(segments[1].a.y, p));
    REQUIRE(segments[1].b.x == 1.0);
    REQUIRE(segments[1].b.y == 0.0);
  } else {
    auto step = std::numbers::pi_v<double> / 2.0 / static_cast<double>(n);
    REQUIRE(segments[0].a.x == 0.0);
    REQUIRE(segments[0].a.y == 1.0);
    for (size_t i = 1; i < n; ++i) {
      auto x = std::pow(std::sin(step * static_cast<double>(i)), 2.0 / d);
      auto y = std::pow(std::cos(step * static_cast<double>(i)), 2.0 / d);
      REQUIRE(apm::details::is_close(segments[i - 1].b.x, x));
      REQUIRE(apm::details::is_close(segments[i - 1].b.y, y));
      REQUIRE(apm::details::is_close(segments[i].a.x, x));
      REQUIRE(apm::details::is_close(segments[i].a.y, y));
    }
    REQUIRE(segments[n - 1].b.x == 1.0);
    REQUIRE(segments[n - 1].b.y == 0.0);
  }
}

TEST_CASE("analytical_model d=1.0", "[apm][analytical_model]") {
  auto model = apm::analytical_model(1.0);
  auto res = std::array{
      2.5e-1, 6.25e-2, 1.5625e-2, 3.90625e-3, 9.765625e-4, 2.441406e-4, 6.103516e-5, 1.525879e-05, 3.814697e-06,
  };
  auto iter1 = model.begin();
  auto iter2 = res.begin();
  for (size_t i = 0, p = 1; i < res.size(); ++i, p <<= 1, ++iter2) {
    for (size_t j = 0; j < p; ++j, ++iter1) {
      REQUIRE(apm::details::is_close(*iter1, *iter2, 1e-5, 0.0));
    }
  }
}

TEST_CASE("analytical_model", "[apm][analytical_model]") {
  auto d = GENERATE(take(100, random(1.0, 10.0)));
  auto p = std::pow(std::sin(std::numbers::pi_v<double> / 4.0), 2.0 / d);
  auto model = apm::analytical_model(d);
  auto iter = model.begin();
  auto rtol = 1e-8;
  auto atol = 0.0;
  for (size_t i = 0; i < 1000; ++i, ++iter) {
    if (i == 0) {
      REQUIRE(apm::details::is_close(*iter, p * p, rtol, atol));
    } else {
      auto res = (1.0 - p) * p / std::pow(4.0, static_cast<int>(std::log2(i + 1)));
      REQUIRE(apm::details::is_close(*iter, res, rtol, atol));
    }
  }
}

TEST_CASE("analytical_model == exact_model == greedy_model", "[apm][analytical_model][exact_model][greedy_model]") {
  auto d = GENERATE(take(100, random(1.0, 10.0)));
  auto model1 = apm::analytical_model(d);
  auto model2 = apm::exact_model(apm::piecewise_segments(2, d), {0.0, 0.0});
  auto model3 = apm::greedy_model(apm::piecewise_segments(2, d), {0.0, 0.0});
  auto iter1 = model1.begin();
  auto iter2 = model2.begin();
  auto iter3 = model3.begin();
  auto rtol = 1e-8;
  auto atol = 0.0;
  for (size_t i = 0; i < 1000; ++i, ++iter1, ++iter2, ++iter3) {
    REQUIRE(apm::details::is_close(*iter1, *iter2, rtol, atol));
    REQUIRE(apm::details::is_close(*iter1, (*iter3).second, rtol, atol));
  }
}

TEST_CASE("exact_model == greedy_model", "[apm][exact_model][greedy_model]") {
  auto d = GENERATE(take(50, random(0.1, 1.0)));
  auto r = GENERATE(chunk(2, take(50, random(-10.0, 0.0))));
  auto model1 = apm::exact_model(apm::piecewise_segments(2, d), {r[0], r[1]});
  auto model2 = apm::greedy_model(apm::piecewise_segments(2, d), {r[0], r[1]});
  auto iter1 = model1.begin();
  auto iter2 = model2.begin();
  auto rtol = 1e-8;
  auto atol = 0.0;
  for (size_t i = 0; i < 1000; ++i, ++iter1, ++iter2) {
    REQUIRE(apm::details::is_close(*iter1, (*iter2).second, rtol, atol));
  }
}
