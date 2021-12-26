#include <benchmark/benchmark.h>

#include "../bsm.h"

static void BM_EuropeanCall(benchmark::State& state) {
    double const K = 100.0;
    var S = 105.0;
    var sigma = 5;
    var variance = sigma*sigma;
    var tau = 30.0 / 365;
    var r = 1.25 / 100;
    var q = 0.0;

    for (auto _: state) {
        var call = european(CP::call, K, S, variance, tau, r, q);
    }
}
BENCHMARK(BM_EuropeanCall);

static void BM_EuropeanCall_Greeks(benchmark::State& state) {
    double const K = 100.0;
    var S = 105.0;
    var sigma = 5;
    var variance = sigma*sigma;
    var tau = 30.0 / 365;
    var r = 1.25 / 100;
    var q = 0.0;

    var call = european(CP::call, K, S, variance, tau, r, q);

    for (auto _: state) {
        auto [delta, vega, theta, rho] = derivativesx(call, wrt(S, sigma, tau, r));
        auto [gamma] = derivativesx(delta, wrt(S));
    }
}
BENCHMARK(BM_EuropeanCall_Greeks);

static void BM_EuropeanCall_IVol(benchmark::State& state) {
    double const K = 100.0;
    var S = 105.0;
    var sigma = 0.30;
    var variance = sigma*sigma;
    var tau = 30.0 / 365;
    var r = 1.25 / 100;
    var q = 0.05;

    var call = european(CP::call, K, S, variance, tau, r, q);
    var iv = implied_volatility(CP::call, K, call, S, tau, r, q);

    for (auto _: state) {
        var iv = implied_volatility(CP::call, K, call, S, tau, r, q);
    }
}
BENCHMARK(BM_EuropeanCall_IVol);

BENCHMARK_MAIN();