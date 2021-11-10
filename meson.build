project('apm', 'cpp', version : '0.1', default_options : ['warning_level=3', 'cpp_std=c++20'])

## Header only library

incdir = include_directories('include')
apm_dep = declare_dependency(include_directories: incdir)
install_subdir('include/apm', install_dir: 'include')
pkg = import('pkgconfig')
pkg.generate(
  name: meson.project_name(),
  filebase: meson.project_name(),
  description: 'Anytime Performance Model library',
  subdirs: meson.project_name(),
)

## Executables

executable('apm', 'apps/apm.cpp', dependencies: [apm_dep], install: true)

## Tests

## Benchmark

bench_deps = [dependency('benchmark'), apm_dep]
bench_srcs = ['bench/apm_benchmark.cpp']
benchmark('apm', executable('apm_benchmark', bench_srcs, dependencies: bench_deps))
