#include <benchmark/benchmark.h>
#include <autodiff/reverse/var.hpp>
#include <autodiff/reverse/var/eigen.hpp>

#include "../analytical.h"
#include "../options.h"
#include "../greeks.h"

using namespace autodiff;

static void EuropeanCall_Pricing(benchmark::State& state) {
    double const K = 100.0;
    var S = 105.0;
    var sigma = 5;
    var tau = 30.0 / 365;
    var r = 1.25 / 100;
    var q = 0.0;

    options::params<var> market_state{ S, sigma, tau, r, q };
    analytical::analyticalmethod<var> bsm_reverse_mode{market_state};

    options::european_call<var> europeanCall{K};

    for (auto _: state) {
        auto callPricing = bsm_reverse_mode.solve(europeanCall);
    }
}
BENCHMARK(EuropeanCall_Pricing);


static void EuropeanCall_Delta(benchmark::State& state) {
    using namespace greeks;
    double const K = 100.0;
    var S = 105.0;
    var sigma = 5;
    var tau = 30.0 / 365;
    var r = 1.25 / 100;
    var q = 0.0;

    options::params<var> market_state{ S, sigma, tau, r, q };
    analytical::analyticalmethod<var> bsm_reverse_mode{market_state};

    options::european_call<var> europeanCall{K};

    auto callPricing = bsm_reverse_mode.solve(europeanCall);
    for (auto _: state) {
        auto delta = greeks::delta(callPricing);
    }
}
BENCHMARK(EuropeanCall_Delta);

static void EuropeanCall_All_Greeks(benchmark::State& state) {
    using namespace greeks;
    double const K = 100.0;
    var S = 105.0;
    var sigma = 5;
    var tau = 30.0 / 365;
    var r = 1.25 / 100;
    var q = 0.0;

    options::params<var> market_state{ S, sigma, tau, r, q };
    analytical::analyticalmethod<var> bsm_reverse_mode{market_state};

    options::european_call<var> europeanCall{K};

    auto callPricing = bsm_reverse_mode.solve(europeanCall);
    for (auto _: state) {
        auto [delta, gamma, vega, theta, rho, psi] = all_greeks(callPricing);
    }
}
BENCHMARK(EuropeanCall_All_Greeks);

static void EuropeanCall_IVol(benchmark::State& state) {
    using namespace greeks;
    double const K = 100.0;
    var S = 105.0;
    var sigma = 0.30;
    var tau = 30.0 / 365;
    var r = 1.25 / 100;
    var q = 0.05;

    options::params<var> market_state{ S, sigma, tau, r, q };
    analytical::analyticalmethod<var> bsm_reverse_mode{market_state};

    options::european_call<var> europeanCall{K};

    auto callPricing = bsm_reverse_mode.solve(europeanCall);
    for (auto _: state) {
        auto iv = bsm_reverse_mode.imply_volatility(europeanCall, val((var)callPricing));
    }
}
BENCHMARK(EuropeanCall_IVol);

BENCHMARK_MAIN();