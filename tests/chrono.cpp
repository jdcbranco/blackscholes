#include <catch2/catch.hpp>

#include "../bsm/bsm.h"

#include <chrono>
#include <iostream>

using namespace bsm;
using namespace bsm::chrono;

TEST_CASE("Datetime") {
    datetime maturity {30d / January / 2022};
    datetime maturity1 = 30d / January / 2023;
    datetime maturity2 {30d / January / 2022, 12h + 0min};
    datetime maturity3 = maturity + 0.5_years;

    auto time1 = time_between(maturity, maturity1);
    CHECK(time1.count() == Approx(365.0 / 365.2425));

    auto time2 = time_between(maturity, maturity2);
    CHECK(time2.count()==Approx(0.5/365.2425));

    auto time3 = time_between(maturity, maturity3);
    CHECK(time3.count()==Approx(0.5));
}


TEST_CASE("Chrono Double literals") {
    auto time = 0.5_years;
    CHECK(time.count()==Approx(0.5));
}

