#include <iapm.hpp>

#include <array>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  if (argc < 6) {
    std::cout << "Error: Missing arguments!\n";
    std::cout << "Usage: " << argv[0] << " n l d r1 r2 [greedy]\n";
    return 1;
  }

  int n = std::stoi(argv[1]);
  int l = std::stoi(argv[2]);
  double d = std::stod(argv[3]);
  std::array<double, 2> r = {std::stod(argv[4]), std::stod(argv[5])};

  int greedy = 0;
  if (argc >= 7) {
    greedy = std::stoi(argv[6]);
  }

  auto hvs = [&]() -> std::vector<double> {
    if (greedy != 0)
      return iapm::iapm_greedy(n, l, d, r, 2);
    else
      return iapm::iapm(n, l, d, r, 2);
  }();

  for (auto const &v : hvs) {
    std::cout << v << "\n";
  }

  return 0;
}
