#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "../options.h"
#include "../analytical.h"
#include "../greeks.h"

using namespace autodiff;

TEST_CASE("Call pricing") {
    double const K = 100.0;
    var S = 105.0;
    var sigma = 5;
    var tau = 30.0 / 365;
    var r = 1.25 / 100;
    var q = 0.0;

    options::params<var> market_state{ S, sigma, tau, r, q };
    analytical::analyticalmethod<var> bsm_reverse_mode{market_state};

    options::european_call<var> europeanCall{K};

    auto [call, greeks] = bsm_reverse_mode.solveForGreeks(europeanCall);
    auto [delta, gamma, vega, theta, rho, psi] = greeks;

    CHECK(call==Approx(56.5136));
    CHECK(delta==Approx(0.773818));
    CHECK(gamma==Approx(0.00199852));
    REQUIRE(vega > 0.0);
}

TEST_CASE("Call Put Parity") {
    double const K = 100.0;
    var S = 105.0;
    var sigma = 5;
    var tau = 30.0 / 365;
    var r = 1.25 / 100;
    var q = 0.0;

    options::params<var> market_state{ S, sigma, tau, r, q };
    analytical::analyticalmethod<var> bsm_reverse_mode{market_state};

    options::european_call<var> europeanCall{K};
    options::european_put<var> europeanPut{K};
    options::forward<var> forwardInstrument{K};

    auto [call, call_greeks] = bsm_reverse_mode.solveForGreeks(europeanCall);
    auto [put, put_greeks] = bsm_reverse_mode.solveForGreeks(europeanPut);

    auto [call_delta, call_gamma, call_vega, call_theta, call_rho, call_psi] = call_greeks;
    auto [put_delta, put_gamma, put_vega, put_theta, put_rho, put_psi] = put_greeks;

    auto fwdPricing = bsm_reverse_mode.solve(forwardInstrument);
    auto [fwd_delta] = greeks::delta(fwdPricing);
    auto fwd = (var)fwdPricing;

    //CHECK((call-put)==Approx(fwd))
    REQUIRE( abs(call-put-fwd) <= 1.0e-5 );
    REQUIRE( abs(call_delta-put_delta-fwd_delta) <= 1.0e-10 );
    REQUIRE( abs(call_gamma-put_gamma) <= 1.0e-10 );
    REQUIRE( abs(call_vega-put_vega) <= 1.0e-10 );

}

