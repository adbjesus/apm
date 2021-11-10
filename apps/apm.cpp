#include <apm/apm.hpp>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  if (argc < 8) {
    std::cout << "Error: Missing arguments!\n";
    std::cout << "Usage: " << argv[0] << " n l d r1 r2 greedy type\n";
    return 1;
  }

  int n = std::stoi(argv[1]);
  int l = std::stoi(argv[2]);
  double d = std::stod(argv[3]);
  double rx = std::stod(argv[4]);
  double ry = std::stod(argv[5]);

  int greedy = std::stoi(argv[6]);
  int type = std::stoi(argv[7]);

  int ulp = 50;

  if (type == 0) {
    auto hvs = [&]() -> std::vector<double> {
      if (greedy == 0)
        return apm::tapm(n, l, d, rx, ry, ulp);
      else
        return apm::tapm_greedy(n, l, d, rx, ry, ulp);
    }();

    for (auto const &v : hvs) {
      std::cout << v << "\n";
    }
  } else {
    auto points = [&]() -> std::vector<apm::point<double>> {
      if (greedy == 0)
        throw("Not yet implemented!");
      else
        return apm::tapm_greedy_points(n, l, d, rx, ry, ulp);
    }();

    for (auto const &p : points) {
      std::cout << p.x << " " << p.y << "\n";
    }
  }

  return 0;
}
