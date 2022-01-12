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

TEST_CASE("European Put Pricing using Binomial Tree (CRR)") {
    auto K = 100.0;
    auto S = 100.0;
    auto sigma = 0.20;
    auto t = system_clock::now();
    auto r = 0.01;
    auto q = 0.05;
    mkt_params mktParams{S, sigma, t, r, q};
    european_put europeanPut{K, t + 0.5_years};
    crr_solver solve{mktParams,2000};
    analytical_solver<autodiff_var> solve_var{mktParams};

    auto crrPricing = solve(europeanPut);
    auto varPricing = solve_var(europeanPut);

    CHECK(crrPricing->price() == Approx(varPricing->price()).margin(0.005));
    CHECK(crrPricing->delta() == Approx(varPricing->delta()).margin(0.00005));
    CHECK(crrPricing->gamma() == Approx(varPricing->gamma()).margin(0.00003));
    CHECK(crrPricing->theta() == Approx(varPricing->theta()).margin(0.005));
    CHECK(crrPricing->vega() == Approx(varPricing->vega()).epsilon(0.005));
    CHECK(crrPricing->rho() == Approx(varPricing->rho()).epsilon(0.005));
    CHECK(crrPricing->psi() == Approx(varPricing->psi()).epsilon(0.005));
}

TEST_CASE("American Call Pricing using Binomial Tree (CRR) without dividends") {
    auto K = 100.0;
    auto S = 100.0;
    auto sigma = 0.20;
    auto t = system_clock::now();
    auto r = 0.01;
    auto q = 0.0;
    mkt_params mktParams{S, sigma, t, r, q};
    american_call americanCall{K, t + 0.5_years};
    european_call europeanCall{K, t + 0.5_years};
    crr_solver solve{mktParams,2000};
    analytical_solver<autodiff_var> solve_var{mktParams};

    auto crrPricing = solve(americanCall);
    auto varPricing = solve_var(europeanCall);

    //These numbers are taken from the algo itself, so not correct, but roughly in line with online calculators
    CHECK(crrPricing->price() == Approx(5.8753208697).margin(0.005));
    CHECK(crrPricing->delta() == Approx(0.5422297484).margin(0.00005));
    CHECK(crrPricing->gamma() == Approx(0.0280617422).margin(0.00003));
    CHECK(crrPricing->theta() == Approx(-6.0958554175).margin(0.005));
    CHECK(crrPricing->vega() == Approx(28.0472071444).epsilon(0.005));
    CHECK(crrPricing->rho() == Approx(24.1767423916).epsilon(0.005));
    CHECK(crrPricing->psi() == Approx(-26.2740657739).epsilon(0.005));

    //Comparision between American and European calls in absence of dividends show that the models produce the same prices, as expected
    CHECK(crrPricing->price() == Approx(varPricing->price()).margin(0.005));
    CHECK(crrPricing->delta() == Approx(varPricing->delta()).margin(0.00005));
    CHECK(crrPricing->gamma() == Approx(varPricing->gamma()).margin(0.00003));
    CHECK(crrPricing->theta() == Approx(varPricing->theta()).margin(0.005));
    CHECK(crrPricing->vega() == Approx(varPricing->vega()).epsilon(0.005));
    CHECK(crrPricing->rho() == Approx(varPricing->rho()).epsilon(0.005));
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

    //These are for q=0.05
    //These numbers are taken from the algo itself, so not correct, but roughly in line with online calculators
    CHECK(crrPricing->price() == Approx(6.5933242703).margin(0.005));
    CHECK(crrPricing->delta() == Approx(-0.5151482623).margin(0.00005));
    CHECK(crrPricing->gamma() == Approx(0.0274551564).margin(0.00003));
    CHECK(crrPricing->theta() == Approx(-7.4856732784).margin(0.005));
    CHECK(crrPricing->vega() == Approx(27.4428949973).epsilon(0.005));
    CHECK(crrPricing->rho() == Approx(-29.049575029).epsilon(0.005));
    CHECK(crrPricing->psi() == Approx(25.7710879181).epsilon(0.005));

}

TEST_CASE("American Put Pricing using Binomial Tree (CRR) using Extra Steps") {
    auto K = 100.0;
    auto S = 100.0;
    auto sigma = 0.20;
    auto t = system_clock::now();
    auto r = 0.01;
    auto q = 0.05;
    mkt_params mktParams{S, sigma, t, r, q};
    american_put americanPut{K, t + 0.5_years};

    crr_solver solve{mktParams,2000,200};
    auto crrPricing = solve(americanPut);

    auto exercise_boundary = crrPricing->exercise_boundary(0.5);
    std::cout << "Exercise boundary = " << exercise_boundary << std::endl;

    //These are for q=0.05
    //These numbers are taken from the algo itself, so not correct, but roughly in line with online calculators
    CHECK(crrPricing->price() == Approx(6.5933242703).margin(0.005));
    CHECK(crrPricing->delta() == Approx(-0.5151482623).margin(0.00005));
    CHECK(crrPricing->gamma() == Approx(0.0274551564).margin(0.00003));
    CHECK(crrPricing->theta() == Approx(-7.4856732784).margin(0.005));
    CHECK(crrPricing->vega() == Approx(27.4428949973).epsilon(0.005));
    CHECK(crrPricing->rho() == Approx(-29.049575029).epsilon(0.005));
    CHECK(crrPricing->psi() == Approx(25.7710879181).epsilon(0.005));

}
TEST_CASE("Test New Superpositioned Binomial Lattice method") {
    auto K = 100.0;
    auto S = 100.0;
    auto sigma = 0.20;
    auto t = system_clock::now();
    auto r = 0.01;
    auto q = 0.05;
    mkt_params mktParams{S, sigma, t, r, q};
    american_put americanPut{K, t + 0.5_years};

    sbl_solver solve{mktParams,200};
    auto sbl_method = solve(americanPut);

    crr_solver solve_crr{mktParams,2000,200};
    auto crr_method = solve_crr(americanPut);

    auto crr_price = crr_method->price();
    auto exercise_boundary = crr_method->exercise_boundary(0.5);
    std::cout << "Price using CRR = " << crr_price << std::endl;
    std::cout << "Exercise boundary using CRR = " << exercise_boundary << std::endl;

}