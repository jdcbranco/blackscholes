#include <iostream>
#include <autodiff/reverse/var.hpp>
#include "bsm.h"

using namespace autodiff;

int main()
{
    //Example of European Call pricing and greek calculation using autodiff
    double const K = 15000.0;
    var S = 15000.0;
    var sigma = 0.21;
    var variance = sigma * sigma;
    var tau = 0.5;
    var r = 0.0;
    var q = 0.25;
    var call = european(CP::call, K, S, variance, tau, r, q);

    auto [delta, vega, theta, rho, psi] = derivativesx(call, wrt(S, sigma, tau, r, q));
    auto [gamma] = derivativesx(delta, wrt(S));
    auto [kappa] = derivativesx(call, wrt(variance));

    std::cout << "call = " << call << std::endl;
    std::cout << "delta = "  << delta  << std::endl;
    std::cout << "gamma = "  << gamma  << std::endl;
    std::cout << "vega = " << vega << std::endl;
    std::cout << "kappa = " << kappa << std::endl;
    std::cout << "theta = " << theta << std::endl;
    std::cout << "rho = " << rho << std::endl;
    std::cout << "psi = " << psi << std::endl;
    std::cout << "\n\n";

    //Example of Implied Dividend calculation using Newton-Raphson method and autodiff
    var fwd = forward(K, S, tau, r, q );
    var impl_div = implied_dividend(K, val(fwd), S, tau, r);

    std::cout << "fwd = " << fwd << std::endl;
    std::cout << std::boolalpha;
    std::cout << "impl div = " << impl_div << " which is " << (abs(impl_div-q) <= 1.0e-4)  << std::endl;

    //Example of Implied Volatility calculation using Newton-Raphson method and autodiff
    var iv = implied_volatility(CP::call, K, val(call), S, tau, r, impl_div);

    std::cout << "impl vol = " << iv << " which is " << (abs(iv-sigma) <= 1.0e-9) << std::endl;


    return 0;
}