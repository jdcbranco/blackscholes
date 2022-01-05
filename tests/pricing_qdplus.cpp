#include <catch2/catch.hpp>

#include "../bsm/bsm.h"

#include <chrono>
#include <string>
#include <iostream>
#include <sstream>

using namespace bsm;
using namespace std::chrono;
using namespace std::chrono_literals;
using namespace bsm::chrono;

TEST_CASE("American Put Pricing using QD+ approximation method matches the Binomial tree" ) {
    auto K = 100.0;
    auto S = 100.0;
    auto sigma = 0.20;
    auto t = system_clock::now();
    auto r = 0.01;
    auto q = 0.05;
    mkt_params mktParams{S, sigma, t, r, q};
    american_put americanPut{K, t + 0.5_years};
    qdplus_solver solve{mktParams};

    auto qdplus_method = solve(americanPut);

    //This number to which we compare to is taken from Binomial method (CRR)
    CHECK(qdplus_method->price() == Approx(6.5933242703).margin(0.005));
}

TEST_CASE("American Put Pricing using QD+ approximation method matches results from the original paper, table 7, page 25 (sigma = 0.20)") {
    auto K = 45.0;
    auto S = 40.0;
    auto sigma = 0.20;
    auto t = system_clock::now();
    auto r = 0.0488;
    auto q = 0.0;
    mkt_params mktParams{S, sigma, t, r, q};
    american_put americanPut{K, t + 0.583_years};
    qdplus_solver solve{mktParams};
    crr_solver solve_crr{mktParams,5000};

    auto qdplus_method = solve(americanPut);
    auto crr_method = solve_crr(americanPut);

    CHECK(qdplus_method->price() == Approx(5.253).margin(0.0005));
    CHECK(qdplus_method->exercise_boundary(0.583) == Approx(37.49).margin(0.005));
}


TEST_CASE("American Put Pricing using QD+ approximation method matches results from the original paper, table 7, page 25 (sigma = 0.20)") {
    auto K = 45.0;
    auto S = 40.0;
    auto sigma = 0.30;
    auto t = system_clock::now();
    auto r = 0.0488;
    auto q = 0.0;
    mkt_params mktParams{S, sigma, t, r, q};
    american_put americanPut{K, t + 0.3333333_years};
    qdplus_solver solve{mktParams};

    auto qdplus_method = solve(americanPut);

    CHECK(qdplus_method->price() == Approx(5.687).margin(0.0005));
    CHECK(qdplus_method->exercise_boundary(0.3333333) == Approx(34.68).margin(0.005));
}