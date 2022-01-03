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

TEST_CASE("Forward pricing") {
    auto K = 100.0;
    auto S = 100.0;
    auto sigma = 0.20;
    auto t = system_clock::now();
    auto r = 0.01;
    auto q = 0.0;
    mkt_params mktParams{S, sigma, t, r, q};
    forward fwd{K, t + 1.0_years};
    analytical_solver solve{mktParams};

    auto fwdPricing = solve(fwd);

    auto tau = 1.0;
    auto expected_price = S * exp(-q*tau) -  K * exp(-r*tau);

    CHECK(fwdPricing->price()==Approx(expected_price));
}

TEST_CASE("Forward pricing using Autodiff Dual") {
    auto K = 100.0;
    auto S = 100.0;
    auto sigma = 0.20;
    auto t = system_clock::now();
    auto r = 0.01;
    auto q = 0.0;
    mkt_params mktParams{S, sigma, t, r, q};
    forward fwd{K, t + 0.75_years};
    analytical_solver<autodiff_dual> solve{mktParams};

    auto fwdPricing = solve(fwd);

    auto tau = 0.75;
    auto expected_price = S * exp(-q*tau) -  K * exp(-r*tau);

    CHECK(fwdPricing->price()==Approx(expected_price));
}

TEST_CASE("European Call Pricing") {
    auto K = 100.0;
    auto S = 100.0;
    auto sigma = 0.20;
    auto t = system_clock::now();
    auto r = 0.01;
    auto q = 0.05;
    mkt_params mktParams{S, sigma, t, r, q};
    european_call europeanCall{K, t + 0.5_years};
    analytical_solver solve{mktParams};

    auto callPricing = solve(europeanCall);

    CHECK(callPricing->price()==Approx(4.62377));
    CHECK(callPricing->delta()==Approx(0.460165));
}

TEST_CASE("European Put Pricing") {
    auto K = 100.0;
    auto S = 100.0;
    auto sigma = 0.20;
    auto t = system_clock::now();
    auto r = 0.02;
    auto q = 0.01;
    mkt_params mktParams{S, sigma, t, r, q};
    european_put europeanPut{K, t + 0.5_years};
    analytical_solver solve{mktParams};

    auto putPricing = solve(europeanPut);

    CHECK(putPricing->price()==Approx(5.3504528757));
    CHECK(putPricing->delta()==Approx(-0.4554818745));
    CHECK(putPricing->gamma()==Approx(0.0279113405));
}

TEST_CASE("European Call Pricing using Autodiff Dual") {
    auto K = 100.0;
    auto S = 100.0;
    auto sigma = 0.20;
    auto t = system_clock::now();
    auto r = 0.01;
    auto q = 0.05;
    mkt_params mktParams{S, sigma, t, r, q};
    european_call europeanCall{K, t + 0.5_years};
    analytical_solver<autodiff_dual> solve{mktParams};

    auto callPricing = solve(europeanCall);

    CHECK(callPricing->price()==Approx(4.62377));
    CHECK(callPricing->delta()==Approx(0.460165));
}

TEST_CASE("European Call Pricing using Autodiff Var") {
    auto K = 100.0;
    auto S = 100.0;
    auto sigma = 0.20;
    auto t = system_clock::now();
    auto r = 0.01;
    auto q = 0.05;
    mkt_params mktParams{S, sigma, t, r, q};
    european_call europeanCall{K, t + 0.5_years};
    analytical_solver<autodiff_var> solve{mktParams};

    auto callPricing = solve(europeanCall);

    CHECK(callPricing->price()==Approx(4.62377));
    CHECK(callPricing->delta()==Approx(0.460165));
}

TEST_CASE("European Call Pricing using different methods must agree") {
    auto K = 100.0;
    auto S = 100.0;
    auto sigma = 0.20;
    auto t = system_clock::now();
    auto r = 0.01;
    auto q = 0.05;
    mkt_params mktParams{S, sigma, t, r, q};
    european_call europeanCall{K, t + 0.5_years};
    analytical_solver<autodiff_off> solve1{mktParams};
    analytical_solver<autodiff_dual> solve2{mktParams};
    analytical_solver<autodiff_var> solve3{mktParams};

    auto callPricing = solve1(europeanCall);
    auto callPricing_dual = solve2(europeanCall);
    auto callPricing_var = solve3(europeanCall);

    //Check Prices against the analytical formula
    CHECK(callPricing_dual->price()==Approx(callPricing->price()));
    CHECK(callPricing_var->price()==Approx(callPricing->price()));

    //Check Prices
    CHECK(callPricing_dual->price()==Approx(callPricing_var->price()));

    //Check Delta
    CHECK(callPricing_dual->delta()==Approx(callPricing_var->delta()));

    //Check Gamma
    CHECK(callPricing_dual->gamma()==Approx(callPricing_var->gamma()));

    //Check Vega
    CHECK(callPricing_dual->vega()==Approx(callPricing_var->vega()));
}

