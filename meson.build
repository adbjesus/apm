project('iapm', ['c', 'cpp'],
  version : '0.1',
  default_options : ['warning_level=3', 'cpp_std=c++20'])

incdir = include_directories('include')

# libiapm = library('iapm', 'src/iapm.cpp',
#   include_directories : incdir,
#   dependencies : deps)

executable('iapm', 'src/main.cpp',
  include_directories : incdir)
