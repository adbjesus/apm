#include <apm/apm.hpp>

#include <array>
#include <iostream>
#include <ranges>
#include <string>

int main(int argc, char **argv) {
  if (argc < 7) {
    std::cout << "Error: Missing arguments!\n";
    std::cout << "Usage: " << argv[0] << " n l d r0 r1 greedy\n";
    return 1;
  }

  auto n = std::stoul(argv[1]);
  auto l = std::stoul(argv[2]);
  auto d = std::stod(argv[3]);
  auto r = apm::point{std::stod(argv[4]), std::stod(argv[5])};
  auto greedy = std::stoi(argv[6]);

  auto segments = apm::piecewise_segments(l, d);

  std::cout << "hv,x,y\n";
  if (greedy == 0) {
    auto model = apm::exact_model(std::move(segments), r);
    for (auto hv : std::ranges::views::take(model, n)) {
      std::cout << hv << ",,\n";
    }
  } else {
    auto model = apm::greedy_model(std::move(segments), r);
    for (auto [point, hv] : std::ranges::views::take(model, n)) {
      std::cout << hv << "," << point.x << "," << point.y << "\n";
    }
  }

  return 0;
}
