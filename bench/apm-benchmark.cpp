#include <benchmark/benchmark.h>

#include <apm/apm.hpp>

#include <ranges>

static void BM_greedy_model(benchmark::State& state) {
  int n = state.range(0);
  int l = state.range(1);
  double d = 2.0 / static_cast<double>(state.range(2));
  auto r = apm::point{0.0, 0.0};
  auto segments = apm::piecewise_segments(l, d);
  for (auto _ : state) {
    auto model = apm::greedy_model(segments, r);
    for (auto [_, hv] : model | std::ranges::views::take(n)) {
      benchmark::DoNotOptimize(hv);
    }
  }
}
BENCHMARK(BM_greedy_model)->ArgsProduct({{10, 100, 1000, 10000, 100000}, {10, 100, 1000}, {4, 1}});

static void BM_exact_model(benchmark::State& state) {
  int n = state.range(0);
  int l = state.range(1);
  double d = 2.0 / static_cast<double>(state.range(2));
  auto r = apm::point{0.0, 0.0};
  auto segments = apm::piecewise_segments(l, d);
  for (auto _ : state) {
    auto model = apm::exact_model(segments, r);
    for (auto hv : model | std::ranges::views::take(n)) {
      benchmark::DoNotOptimize(hv);
    }
  }
}
BENCHMARK(BM_exact_model)->ArgsProduct({{10, 100, 1000, 10000, 100000}, {10, 100, 1000}, {4, 1}});

BENCHMARK_MAIN();
