cmake_minimum_required(VERSION 3.18)

if(NOT DEFINED PROJECT_NAME)
  set(APM_MASTER_PROJECT ON)
else()
  set(APM_MASTER_PROJECT OFF)
endif()

project(
  APM
  VERSION 0.1.1
  DESCRIPTION "Anytime Performance Model for bi-objective optimization algorithms."
  LANGUAGES CXX
)

# Project options
option(APM_BUILD_APPS "Build the apps" ${APM_MASTER_PROJECT})
option(APM_BUILD_DOCS "Build the documentation" ${APM_MASTER_PROJECT})
option(APM_BUILD_TESTS "Build tests" OFF)
option(APM_BUILD_BENCHMARKS "Build benchmarks" OFF)
option(APM_ENABLE_WARNINGS "Enable warnings" OFF)
option(APM_ENABLE_WERROR "Enable -Werror or /WX" OFF)

# Warning flags
set(APM_CXX_WARN_FLAGS "")

if (APM_ENABLE_WARNINGS)
  if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    list(APPEND APM_CXX_WARN_FLAGS /W4)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    list(APPEND APM_CXX_WARN_FLAGS
      -Wall
      -Wextra
      -Wpedantic
      -pedantic-errors
      -Wold-style-cast
      -Wundef
      -Wredundant-decls
      -Wcast-qual
      -Wmissing-include-dirs
      -Wdisabled-optimization
      -Wconversion
      -Wdouble-promotion
      -Wshadow
      -Wnull-dereference
    )
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
      list(APPEND APM_CXX_WARN_FLAGS
        -Wvector-operation-performance
        -Wduplicated-cond
      )
    endif()
  else()
    message(WARNING "Unknown compiler (not GNU, Clang, or MSVC), not setting any warning flags.")
    set(APM_CXX_WARN_FLAGS "")
  endif()
endif()

if(APM_ENABLE_WERROR)
  list(APPEND APM_CXX_WARN_FLAGS
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Werror>
    $<$<CXX_COMPILER_ID:MSVC>:/WX>
  )
endif()

# Extra stuff if this is the master project
if(APM_MASTER_PROJECT)
  # Turn off C++ extensions
  set(CMAKE_CXX_EXTENSIONS OFF)
  # Create compile commands json for clang tools
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
  # Include CTest
  include(CTest)
endif()

# Install dirs for install targets
include(GNUInstallDirs)

# Define library target
add_library(apm INTERFACE)
target_compile_features(apm INTERFACE cxx_std_20)
target_include_directories(apm INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# The library is an install target
install(TARGETS apm)
install(DIRECTORY "${PROJECT_SOURCE_DIR}/include/apm" DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# PkgConfig file
configure_file(
  "${PROJECT_SOURCE_DIR}/cmake/apm.pc.in"
  "${PROJECT_BINARY_DIR}/apm.pc"
  @ONLY
)
install(
  FILES "${PROJECT_BINARY_DIR}/apm.pc"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig"
)

# Export
install(
  TARGETS apm
  EXPORT apm-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
  EXPORT apm-targets
  NAMESPACE apm::
  FILE apm-targets.cmake
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/apm"
)

include(CMakePackageConfigHelpers)
configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/apm-config.cmake.in"
  "${PROJECT_BINARY_DIR}/apm-config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/apm"
)
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/apm-config-version.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY SameMinorVersion
)

install(
  FILES "${PROJECT_BINARY_DIR}/apm-config.cmake"
        "${PROJECT_BINARY_DIR}/apm-config-version.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/apm"
)

if(APM_BUILD_APPS)
  add_executable(apm-bin apps/apm.cpp)
  target_link_libraries(apm-bin apm)
  target_compile_options(apm-bin PRIVATE ${APM_CXX_WARN_FLAGS})
  set_target_properties(apm-bin PROPERTIES OUTPUT_NAME apm)
  install(TARGETS apm-bin)
endif()

# Documentation
if(APM_BUILD_DOCS)
  find_package(Doxygen REQUIRED)
  set(DOXYGEN_COLLABORATION_GRAPH NO)
  set(DOXYGEN_QUIET YES)
  set(DOXYGEN_WARN_IF_UNDOCUMENTED NO)
  doxygen_add_docs(docs "${PROJECT_SOURCE_DIR}/include" ALL)
  install(
    DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/html"
    DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/doc/apm"
  )
endif()

if(APM_BUILD_TESTS)
  find_package(Catch2 REQUIRED)
  add_executable(apm-test tests/main.cpp tests/apm/apm.cpp)
  target_link_libraries(apm-test Catch2::Catch2 apm)
  target_compile_options(apm-test PRIVATE ${APM_CXX_WARN_FLAGS})
  include(Catch)
  catch_discover_tests(apm-test)
endif()

# Benchmarks
if(APM_BUILD_BENCHMARKS)
  find_package(benchmark REQUIRED)
  add_executable(apm-benchmark bench/apm-benchmark.cpp)
  target_compile_options(apm-benchmark PRIVATE ${APM_CXX_WARN_FLAGS})
  target_link_libraries(apm-benchmark apm benchmark::benchmark)
endif()
