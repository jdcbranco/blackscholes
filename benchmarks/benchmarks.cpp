#include <benchmark/benchmark.h>

#include "bsm.h"

using namespace bsm;
using namespace bsm::chrono;

static void Benchmark_EC_Baseline_Price(benchmark::State& state) {
    auto K = 100.0;
    auto S = 100.0;
    auto sigma = 0.20;
    auto t = datetime::now();
    auto r = 0.01;
    auto q = 0.05;
    mkt_params mktParams{S, sigma, t, r, q};
    european_call europeanCall{K, t + 0.5_years};
    analytical_solver<autodiff_off> solve{mktParams};

    for (auto _: state) {
        auto callPricing = solve(europeanCall);
        callPricing->price();
    }
}
BENCHMARK(Benchmark_EC_Baseline_Price);

static void Benchmark_EC_Baseline_Delta(benchmark::State& state) {
    auto K = 100.0;
    auto S = 100.0;
    auto sigma = 0.20;
    auto t = datetime::now();
    auto r = 0.01;
    auto q = 0.05;
    mkt_params mktParams{S, sigma, t, r, q};
    european_call europeanCall{K, t + 0.5_years};
    analytical_solver<autodiff_off> solve{mktParams};

    auto callPricing = solve(europeanCall);

    for (auto _: state) {
        callPricing->delta();
    }
}
BENCHMARK(Benchmark_EC_Baseline_Delta);

static void Benchmark_EC_Dual_Price(benchmark::State& state) {
    auto K = 100.0;
    auto S = 100.0;
    auto sigma = 0.20;
    auto t = datetime::now();
    auto r = 0.01;
    auto q = 0.05;
    mkt_params mktParams{S, sigma, t, r, q};
    european_call europeanCall{K, t + 0.5_years};
    analytical_solver<autodiff_dual> solve{mktParams};

    for (auto _: state) {
        auto callPricing = solve(europeanCall);
        callPricing->price();
    }
}
BENCHMARK(Benchmark_EC_Dual_Price);

//Autodiff Reverse mode

static void Benchmark_EC_Var_Price(benchmark::State& state) {
    auto K = 100.0;
    auto S = 100.0;
    auto sigma = 0.20;
    auto t = datetime::now();
    auto r = 0.01;
    auto q = 0.05;
    mkt_params mktParams{S, sigma, t, r, q};
    european_call europeanCall{K, t + 0.5_years};
    analytical_solver<autodiff_var> solve{mktParams};

    for (auto _: state) {
        auto callPricing = solve(europeanCall);
        callPricing->price();
    }
}
BENCHMARK(Benchmark_EC_Var_Price);

//Binomial method

static void Benchmark_EC_CRR_Price(benchmark::State& state) {
    auto K = 100.0;
    auto S = 100.0;
    auto sigma = 0.20;
    auto t = datetime::now();
    auto r = 0.01;
    auto q = 0.05;
    mkt_params mktParams{S, sigma, t, r, q};
    european_call europeanCall{K, t + 0.5_years};
    crr_solver<autodiff_off> solve{mktParams,400};

    for (auto _: state) {
        auto pricing = solve(europeanCall);
        pricing->price();
        pricing->delta();
        pricing->gamma();
        pricing->vega();
        pricing->rho();
        pricing->theta();
        pricing->psi();
    }
}
BENCHMARK(Benchmark_EC_CRR_Price);

static void Benchmark_AP_CRR_Price(benchmark::State& state) {
    auto K = 100.0;
    auto S = 100.0;
    auto sigma = 0.20;
    auto t = datetime::now();
    auto r = 0.01;
    auto q = 0.05;
    mkt_params mktParams{S, sigma, t, r, q};
    american_put americanPut{K, t + 0.5_years};
    crr_solver<autodiff_off> solve{mktParams,400};

    for (auto _: state) {
        auto pricing = solve(americanPut);
        pricing->price();
        pricing->delta();
        pricing->gamma();
        pricing->vega();
        pricing->rho();
        pricing->theta();
        pricing->psi();
    }
}
BENCHMARK(Benchmark_AP_CRR_Price);

BENCHMARK_MAIN();