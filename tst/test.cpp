#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "../bsm.h"

var kappa(double K, var S, var sigma, var tau, var r) {
    auto const d1 = (log(S / K) + (r + sigma * sigma / 2) * tau) / (sigma * sqrt(tau));
    return 0.5 * S * sqrt(tau) * one_div_root_two_pi * exp(-0.5 * d1*d1) / sigma;
}

TEST_CASE("Call pricing") {
    double const K = 100.0;
    var S = 105.0;
    var sigma = 5;
    var variance = sigma*sigma;
    var tau = 30.0 / 365;
    var r = 1.25 / 100;
    var q = 0.0;
    var call = european(CP::call, K, S, variance, tau, r, q);

    auto [delta, vega, theta, rho] = derivativesx(call, wrt(S, sigma, tau, r));
    auto [gamma] = derivativesx(delta, wrt(S));

    REQUIRE( abs(call-56.5136) <= 1.0e-5);
    REQUIRE( abs(delta-0.773818) <= 1.0e-5);
    REQUIRE( abs(gamma-0.00199852) <= 1.0e-8);
}

TEST_CASE("Call Put Parity") {
    double const K = 100.0;
    var S = 105.0;
    var sigma = 5;
    var variance = sigma*sigma;
    var tau = 30.0 / 365;
    var r = 1.25 / 100;
    var q = 0.0;
    var call = european(CP::call, K, S, variance, tau, r, q);
    var put = european(CP::put, K, S, variance, tau, r, q);

    auto [call_delta, call_vega, call_theta, call_rho] = derivativesx(call, wrt(S, sigma, tau, r));
    auto [call_gamma] = derivativesx(call_delta, wrt(S));

    auto [put_delta, put_vega, put_theta, put_rho] = derivativesx(put, wrt(S, sigma, tau, r));
    auto [put_gamma] = derivativesx(put_delta, wrt(S));

    var fwd = forward(K, S, tau, r, q);
    auto [fwd_delta] = derivativesx(fwd, wrt(S));

    REQUIRE( abs(call-put-fwd) <= 1.0e-5 );
    REQUIRE( abs(call_delta-put_delta-fwd_delta) <= 1.0e-10 );
    REQUIRE( abs(call_gamma-put_gamma) <= 1.0e-10 );
    REQUIRE( abs(call_vega-put_vega) <= 1.0e-10 );

}

TEST_CASE("Test kappa") {
    double const K = 15000.0;
    var S = 15000.0;
    var sigma = 0.21;
    var variance = sigma * sigma;
    var tau = 0.5;
    var r = 0.0;
    var q = 0.0;
    var call = european(CP::call, K, S, variance, tau, r, q);

    auto [autodiff_kappa] = derivativesx(call, wrt(variance));
    auto calculated_kappa = kappa(K, S, sigma, tau, r);

    REQUIRE(abs(autodiff_kappa - calculated_kappa) <= 1.0e-10 );

}