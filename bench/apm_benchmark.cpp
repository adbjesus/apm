#include <benchmark/benchmark.h>

#include <apm/apm.hpp>
#include <iostream>

static void BM_tapm(benchmark::State& state) {
  for (auto _ : state) {
    benchmark::DoNotOptimize(apm::tapm<double>(state.range(0), state.range(1), 2.0 / double(state.range(2)), 0, 0));
  }
}
// Register the function as a benchmark
BENCHMARK(BM_tapm)->ArgsProduct({{10, 100, 1000, 10000, 100000}, {10, 100, 1000}, {4, 1}});

// Define another benchmark
static void BM_tapm_greedy(benchmark::State& state) {
  for (auto _ : state) {
    benchmark::DoNotOptimize(
        apm::tapm_greedy<double>(state.range(0), state.range(1), 2.0 / double(state.range(2)), 0, 0));
  }
}
BENCHMARK(BM_tapm_greedy)->ArgsProduct({{10, 100, 1000, 10000, 100000}, {10, 100, 1000}, {4, 1}});

BENCHMARK_MAIN();
