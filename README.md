# Anytime Performance Model (APM)

A C++ library and binary for computing a models of anytime
performance. Currently, it supports a theoretical model of anytime
performance for algorithms that find at each iteration a solution to a
bi-objective problem, which was first described in:

A. D. Jesus, L. Paquete, A. Liefooghe. **A model of anytime algorithm performance
for bi-objective optimization.** Journal of Global Optimization 78, 329–350,
2021. [DOI](https://doi.org/10.1007/s10898-020-00909-9)

Note: This library is still under active development and breaking
changes are likelly to occur between minor versions. Nonetheless, the
code is tested and all algorithms and data structures implemented are
expected to be correct.

## Dependencies

- [Catch2](https://github.com/catchorg/Catch2) for tests.
- [Doxygen](https://www.doxygen.nl/index.html) for generating documentation.
- [benchmark](https://github.com/google/benchmark) for benchmarks.
