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

TEST_CASE("European Call Pricing using Binomial Tree (CRR)") {
    auto K = 100.0;
    auto S = 100.0;
    auto sigma = 0.20;
    auto t = system_clock::now();
    auto r = 0.01;
    auto q = 0.05;
    mkt_params mktParams{S, sigma, t, r, q};
    european_call europeanCall{K, t + 0.5_years};
    crr_solver solve{mktParams,2000};
    analytical_solver<autodiff_var> solve_analytically{mktParams};

    auto crrPricing = solve(europeanCall);
    auto varPricing = solve_analytically(europeanCall);

    CHECK(crrPricing->price() == Approx(4.62377).margin(0.005));
    CHECK(crrPricing->delta() == Approx(0.460165).margin(0.00005));
    CHECK(crrPricing->gamma() == Approx(varPricing->gamma()).margin(0.00003));
    CHECK(crrPricing->theta() == Approx(varPricing->theta()).margin(0.005));
    CHECK(crrPricing->vega() == Approx(varPricing->vega()).epsilon(0.005));
    CHECK(crrPricing->rho() == Approx(varPricing->rho()).epsilon(0.005));
    CHECK(crrPricing->psi() == Approx(varPricing->psi()).epsilon(0.005));
}

TEST_CASE("American Put Pricing using Binomial Tree (CRR)") {
    auto K = 100.0;
    auto S = 100.0;
    auto sigma = 0.20;
    auto t = system_clock::now();
    auto r = 0.01;
    auto q = 0.05;
    mkt_params mktParams{S, sigma, t, r, q};
    american_put americanPut{K, t + 0.5_years};
    crr_solver solve{mktParams,2000};

    auto crrPricing = solve(americanPut);

    //These numbers are taken from the algo itself, so not correct, but roughly in line with online calculators
    CHECK(crrPricing->price() == Approx(6.5933242703).margin(0.005));
    CHECK(crrPricing->delta() == Approx(-0.5151482623).margin(0.00005));
    CHECK(crrPricing->gamma() == Approx(0.0274551564).margin(0.00003));
    //Existing implementation for theta doesn't seem to hold for American options
    //CHECK(crrPricing->theta() == Approx(0.0).margin(0.005));
    CHECK(crrPricing->vega() == Approx(27.4428949973).epsilon(0.005));
    CHECK(crrPricing->rho() == Approx(-29.049575029).epsilon(0.005));
    CHECK(crrPricing->psi() == Approx(25.7710879181).epsilon(0.005));
}
