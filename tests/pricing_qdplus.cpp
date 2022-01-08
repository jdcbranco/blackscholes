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

TEST_CASE("American Call Pricing using QD+ approximation method matches the Binomial Tree (without dividends)") {
    auto K = 100.0L;
    auto S = 100.0L;
    auto sigma = 0.20L;
    auto t = system_clock::now();
    auto r = 0.01L;
    auto q = 0.0L;
    mkt_params<long double> mktParams{S, sigma, t, r, q};
    american_call americanCall{K, t + 0.5_years};

    qdplus_solver solve{mktParams};
    auto qdplus_method = solve(americanCall);

    crr_solver solve_crr{mktParams,2000};
    auto crr_method = solve_crr(americanCall);

    CHECK(qdplus_method->price() == Approx(crr_method->price()).margin(0.005));
    CHECK(qdplus_method->delta() == Approx(crr_method->delta()).margin(0.0005));
    CHECK(qdplus_method->gamma() == Approx(crr_method->gamma()).margin(0.00003));
    CHECK(qdplus_method->theta() == Approx(crr_method->theta()).epsilon(0.005));
    CHECK(qdplus_method->vega() == Approx(crr_method->vega()).epsilon(0.005));
    CHECK(qdplus_method->rho() == Approx(crr_method->rho()).epsilon(0.005));
    CHECK(qdplus_method->psi() == Approx(crr_method->psi()).epsilon(0.035));
}

//Note that QD+ diverges a lot when diviends are higher
TEST_CASE("American Call Pricing using QD+ approximation method matches the Binomial Tree (with small dividends)") {
    auto K = 100.0L;
    auto S = 100.0L;
    auto sigma = 0.20L;
    auto t = system_clock::now();
    auto r = 0.01L;
    auto q = 0.005L;
    mkt_params<long double> mktParams{S, sigma, t, r, q};
    american_call americanCall{K, t + 0.5_years};
    qdplus_solver solve{mktParams};
    auto qdplus_method = solve(americanCall);

    crr_solver solve_crr{mktParams,2000};
    auto crr_method = solve_crr(americanCall);

    //This number to which we compare to is taken from Binomial method (CRR)
    CHECK(qdplus_method->price() == Approx(crr_method->price()).margin(0.005));
    CHECK(qdplus_method->delta() == Approx(crr_method->delta()).margin(0.0005));
    CHECK(qdplus_method->gamma() == Approx(crr_method->gamma()).margin(0.00003));
    CHECK(qdplus_method->theta() == Approx(crr_method->theta()).epsilon(0.005));
    CHECK(qdplus_method->vega() == Approx(crr_method->vega()).epsilon(0.005));
    CHECK(qdplus_method->rho() == Approx(crr_method->rho()).epsilon(0.005));
    CHECK(qdplus_method->psi() == Approx(crr_method->psi()).epsilon(0.005));

}

TEST_CASE("American Put Pricing using QD+ approximation method matches the Binomial tree" ) {
    auto K = 100.0L;
    auto S = 100.0L;
    auto sigma = 0.20L;
    auto t = system_clock::now();
    auto r = 0.01L;
    auto q = 0.05L;
    mkt_params<long double> mktParams{S, sigma, t, r, q};
    american_put americanPut{K, t + 0.5_years};
    qdplus_solver solve{mktParams};

    auto qdplus_method = solve(americanPut);

    //This number to which we compare to is taken from Binomial method (CRR)
    CHECK(qdplus_method->price() == Approx(6.5933242703).margin(0.005));
    CHECK(qdplus_method->delta() == Approx(-0.5151482623).margin(0.00005));
    CHECK(qdplus_method->gamma() == Approx(0.0274551564).margin(0.00003));
    CHECK(qdplus_method->theta() == Approx(-7.4856732784).margin(0.005));
    CHECK(qdplus_method->vega() == Approx(27.4428949973).epsilon(0.005));
    CHECK(qdplus_method->rho() == Approx(-29.049575029).epsilon(0.005));
    CHECK(qdplus_method->psi() == Approx(25.7710879181).epsilon(0.03));
}

TEST_CASE("American Put Pricing using QD+ approximation method matches results from the original paper, table 7, page 25 (sigma = 0.20)") {
    auto K = 45.0;
    auto S = 40.0;
    auto sigma = 0.20;
    auto t = system_clock::now();
    auto r = 0.0488;
    auto q = 0.0;
    mkt_params<long double> mktParams{S, sigma, t, r, q};
    american_put americanPut{K, t + 0.583_years};
    qdplus_solver solve{mktParams};
    crr_solver solve_crr{mktParams,5000};

    auto qdplus_method = solve(americanPut);
    auto crr_method = solve_crr(americanPut);

    CHECK(qdplus_method->price() == Approx(5.253).margin(0.0005));
    CHECK(qdplus_method->exercise_boundary(0.583) == Approx(37.49).margin(0.005));
}


TEST_CASE("American Put Pricing using QD+ approximation method matches results from the original paper, table 7, page 25 (sigma = 0.30)") {
    auto K = 45.0;
    auto S = 40.0;
    auto sigma = 0.30;
    auto t = system_clock::now();
    auto r = 0.0488;
    auto q = 0.0;
    mkt_params<long double> mktParams{S, sigma, t, r, q};
    american_put americanPut{K, t + 0.3333333_years};
    qdplus_solver solve{mktParams};

    auto qdplus_method = solve(americanPut);

    CHECK(qdplus_method->price() == Approx(5.687).margin(0.0005));
    CHECK(qdplus_method->exercise_boundary(0.3333333) == Approx(34.68).margin(0.005));
}