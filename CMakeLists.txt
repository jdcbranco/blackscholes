cmake_minimum_required(VERSION 3.20)
project(bsm)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_INSTALL_PREFIX  ${CMAKE_CURRENT_SOURCE_DIR} CACHE PATH "install folder" FORCE)

find_package(Catch2 2.13.7 REQUIRED)
find_package(autodiff REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(benchmark REQUIRED)
find_package(Taskflow 3.3.0 REQUIRED)
find_package(Threads REQUIRED)

#bsm library
add_library(bsm STATIC main.cpp random.cpp random.h bsm/bsm.h bsm/instruments.cpp bsm/instruments.h bsm/solver.h bsm/chrono.h bsm/chrono.cpp bsm/solver_analytical.cpp bsm/solver_analytical_autodiff_dual.cpp bsm/solver_analytical_autodiff_var.cpp bsm/bintree.h bsm/solver_crr.cpp bsm/solver_crr_internals.h bsm/solver_fastamerican.cpp bsm/solver_qdplus.cpp bsm/solver_analytical_internals.h bsm/solver_american_internals.h)
target_include_directories(bsm PRIVATE eigen3 bsm)
target_link_libraries(bsm autodiff::autodiff Eigen3::Eigen)

#Main app
add_executable(main main.cpp)
target_include_directories(main PRIVATE bsm)
target_link_libraries(main bsm Threads::Threads)

#Unit tests
add_executable(tests tests/common.cpp random.cpp random.h bsm/bsm.h tests/instruments.cpp tests/pricing_analytical.cpp tests/chrono.cpp tests/pricing_crr.cpp tests/pricing_qdplus.cpp)
target_include_directories(tests PRIVATE eigen3 bsm)
target_link_libraries(tests bsm autodiff::autodiff Eigen3::Eigen )

#Benchmarks
add_executable(benchmarks benchmarks/benchmarks.cpp random.cpp random.h bsm/bsm.h)
target_include_directories(benchmarks PRIVATE eigen3 bsm)
target_link_libraries(benchmarks bsm autodiff::autodiff Eigen3::Eigen benchmark::benchmark)

install(TARGETS bsm DESTINATION bin)