#include <catch2/catch.hpp>

#include "../bsm/bsm.h"

#include <chrono>

using namespace bsm;
using namespace bsm::chrono;
using namespace std::chrono_literals;
using std::chrono::January;
using std::chrono::system_clock;
using std::chrono::sys_days;

TEST_CASE("Forward") {
    auto now = system_clock::now();
    auto oneYearFromNow = now + 1.0_years;

    forward fwd{100.0, oneYearFromNow};
    CHECK(fwd.K==Approx(100.0));
    CHECK(fwd.maturity==oneYearFromNow);

    auto value = fwd.payoff(100.0);
    CHECK(value == Approx(0.0));
}

TEST_CASE("European Call") {
    auto strike = 100.0;
    auto maturity = 30d/January/2022;

    european_call europeanCall{strike, maturity};
    CHECK(europeanCall.K==Approx(strike));
    //CHECK(europeanCall.maturity==sys_days(maturity));
    CHECK(europeanCall.maturity==maturity);

    auto price_at_maturity_itm = 101.0;
    auto price_at_maturity_atm = 100.0;
    auto price_at_maturity_otm = 99.0;

    auto value_at_maturity_itm = europeanCall.payoff(price_at_maturity_itm);
    auto value_at_maturity_atm = europeanCall.payoff(price_at_maturity_atm);
    auto value_at_maturity_otm = europeanCall.payoff(price_at_maturity_otm);

    CHECK(value_at_maturity_itm==Approx(price_at_maturity_itm-europeanCall.K));
    CHECK(value_at_maturity_atm==Approx(0.0));
    CHECK(value_at_maturity_otm==Approx(0.0));
}

TEST_CASE("European Put") {
    auto strike = 100.0;
    auto maturity = 30d/January/2023;

    european_put europeanPut{strike, maturity};

    CHECK(europeanPut.K==Approx(strike));
    CHECK(europeanPut.maturity==maturity);

    auto price_at_maturity_itm = 99.0;
    auto price_at_maturity_atm = 100.0;
    auto price_at_maturity_otm = 101.0;

    auto value_at_maturity_itm = europeanPut.payoff(price_at_maturity_itm);
    auto value_at_maturity_atm = europeanPut.payoff(price_at_maturity_atm);
    auto value_at_maturity_otm = europeanPut.payoff(price_at_maturity_otm);

    CHECK(value_at_maturity_itm==Approx(europeanPut.K-price_at_maturity_itm));
    CHECK(value_at_maturity_atm==Approx(0.0));
    CHECK(value_at_maturity_otm==Approx(0.0));

}

TEST_CASE("American Call") {
    auto strike = 100.0;
    auto maturity = 30d/January/2022;

    american_call americanCall{strike, maturity};
    CHECK(americanCall.K==Approx(strike));
    CHECK(americanCall.maturity==maturity);

    auto price_at_maturity_itm = 101.0;
    auto price_at_maturity_atm = 100.0;
    auto price_at_maturity_otm = 99.0;

    auto value_at_maturity_itm = americanCall.payoff(price_at_maturity_itm);
    auto value_at_maturity_atm = americanCall.payoff(price_at_maturity_atm);
    auto value_at_maturity_otm = americanCall.payoff(price_at_maturity_otm);

    CHECK(value_at_maturity_itm==Approx(price_at_maturity_itm-americanCall.K));
    CHECK(value_at_maturity_atm==Approx(0.0));
    CHECK(value_at_maturity_otm==Approx(0.0));
}

TEST_CASE("American Put") {
    auto strike = 100.0;
    auto maturity = 30d/January/2022;

    american_put americanPut{strike, maturity};
    CHECK(americanPut.K==Approx(strike));
    CHECK(americanPut.maturity==maturity);

    auto price_at_maturity_itm = 98.0;
    auto price_at_maturity_atm = 100.0;
    auto price_at_maturity_otm = 102.0;

    auto value_at_maturity_itm = americanPut.payoff(price_at_maturity_itm);
    auto value_at_maturity_atm = americanPut.payoff(price_at_maturity_atm);
    auto value_at_maturity_otm = americanPut.payoff(price_at_maturity_otm);

    CHECK(value_at_maturity_itm==Approx(americanPut.K-price_at_maturity_itm));
    CHECK(value_at_maturity_atm==Approx(0.0));
    CHECK(value_at_maturity_otm==Approx(0.0));
}