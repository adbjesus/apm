#include <apm/theoretical_model.hpp>

#include <array>
#include <iostream>
#include <ranges>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  if (argc < 7) {
    std::cout << "Error: Missing arguments!\n";
    std::cout << "Usage: " << argv[0] << " n l d r0 r1 greedy\n";
    return 1;
  }

  int n = std::stoi(argv[1]);
  int l = std::stoi(argv[2]);
  double d = std::stod(argv[3]);
  double r0 = std::stod(argv[4]);
  double r1 = std::stod(argv[5]);
  int greedy = std::stoi(argv[6]);

  // std::vector<std::pair<std::optional<apm::point<double>>, double>> result;
  std::vector<double> result;
  result.reserve(n);
  auto segments = apm::piecewise_segments(l, d);

  if (greedy == 0) {
    auto model = apm::exact_model(std::move(segments), {r0, r1});
    for (auto hv : model | std::ranges::views::take(n)) {
      result.emplace_back(hv);
    }
  } else {
    auto model = apm::greedy_model(std::move(segments), {r0, r1});
    for (auto [_, hv] : model | std::ranges::views::take(n)) {
      result.emplace_back(hv);
    }
    // apm::greedy_model(std::back_inserter(hvs), n, l, d, rx, ry);
  }

  for (auto const &hv : result) {
    std::cout << hv << "\n";
  }

  return 0;
}
