cmake_minimum_required(VERSION 3.20)
project(blackscholes)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_INSTALL_PREFIX  ${CMAKE_CURRENT_SOURCE_DIR} CACHE PATH "install folder" FORCE)

find_package(Catch2 2.13.7 REQUIRED)
find_package(autodiff REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(benchmark REQUIRED)

#Main app
add_executable(blackscholes main.cpp crr.cpp crr.h lattice.h analytical.h greeks.h options.cpp options.h random.cpp random.h lsmc.cpp lsmc.h)
target_include_directories(blackscholes PRIVATE eigen3)
target_link_libraries(blackscholes autodiff::autodiff Eigen3::Eigen)


#Unit tests
add_executable(unit_tests tst/test.cpp crr.cpp crr.h lattice.h analytical.h greeks.h options.cpp options.h random.cpp random.h lsmc.cpp lsmc.h)
target_include_directories(unit_tests PRIVATE eigen3)
target_link_libraries(unit_tests autodiff::autodiff Eigen3::Eigen)

#Benchmarks
add_executable(benchmarks benchmarks/benchmarks.cpp crr.cpp crr.h lattice.h analytical.h greeks.h options.cpp options.h random.cpp random.h lsmc.cpp lsmc.h)
target_include_directories(benchmarks PRIVATE eigen3)
target_link_libraries(benchmarks autodiff::autodiff Eigen3::Eigen benchmark::benchmark)

install(TARGETS blackscholes DESTINATION bin)