cmake_minimum_required(VERSION 3.20)
project(blackscholes)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_INSTALL_PREFIX  ${CMAKE_CURRENT_SOURCE_DIR} CACHE PATH "install folder" FORCE)

find_package(Catch2 2.13.7 REQUIRED)
find_package(autodiff REQUIRED)
find_package(benchmark REQUIRED)


#Main app
add_executable(blackscholes main.cpp bsm.cpp bsm.h)
target_link_libraries(blackscholes autodiff::autodiff)

#Unit tests
add_executable(unit_tests tst/test.cpp bsm.cpp bsm.h)
target_link_libraries(unit_tests autodiff::autodiff)

#Benchmarks
add_executable(benchmarks benchmarks/benchmarks.cpp bsm.cpp bsm.h)
target_link_libraries(benchmarks autodiff::autodiff benchmark::benchmark)

install(TARGETS blackscholes DESTINATION bin)