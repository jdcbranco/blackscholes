#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "../bsm/bsm.h"

#include <chrono>

TEST_CASE("Mkt Params") {
    using namespace std::chrono;
    using namespace std::chrono_literals;
    auto S = 100.0;
    auto sigma = 0.20;
    auto t = system_clock::now();
    auto r = 0.01;
    auto q = 0.0;
    bsm::mkt_params mktParams{S, sigma, t, r, q};
    CHECK(mktParams.S == S);
    CHECK(mktParams.sigma == sigma);
    CHECK(mktParams.t == t);
    CHECK(mktParams.r == r);
    CHECK(mktParams.q == q);
}


