#ifndef BLACKSCHOLES_GREEKS_H
#define BLACKSCHOLES_GREEKS_H

#include <autodiff/reverse/var.hpp>
#include <autodiff/reverse/var/eigen.hpp>
#include "options.h"

namespace greeks {
    using namespace autodiff;
    using namespace std;

    template<typename T>
    T delta(options::pricing<T> pricing);
    template<typename T>
    T gamma(options::pricing<T> pricing);
    template<typename T>
    T vega(options::pricing<T> pricing);
    template<typename T>
    T theta(options::pricing<T> pricing);
    template<typename T>
    T rho(options::pricing<T> pricing);
    template<typename T>
    T psi(options::pricing<T> pricing);
    template<typename T>
    T kappa(options::pricing<T> pricing);


    auto delta(options::pricing<var> & pricing) {
        return derivativesx((var&)pricing, wrt(pricing.S()))[0];
    }

    auto gamma(options::pricing<var> & pricing) {
        return derivativesx(delta(pricing), wrt(pricing.S()))[0];
    }

    auto vega(options::pricing<var> & pricing) {
        return derivativesx((var)pricing, wrt(pricing.sigma()))[0];
    }

    auto theta(options::pricing<var> & pricing) {
        return derivativesx((var)pricing, wrt(pricing.tau()))[0];
    }

    auto rho(options::pricing<var> & pricing) {
        return derivativesx((var)pricing, wrt(pricing.r()))[0];
    }

    auto psi(options::pricing<var> & pricing) {
        return derivativesx((var)pricing, wrt(pricing.q()))[0];
    }

    auto all_greeks(options::pricing<var> & pricing) {
        auto& [S, sigma, tau, r, q] = pricing.variables();
        auto [delta, vega, theta, rho, psi] =  derivativesx((var)pricing, wrt(S, sigma, tau, r, q));
        auto gamma = derivativesx(delta, wrt(S))[0];
        array<var, 6> all_greeks{delta, gamma, vega, theta, rho, psi};
        return all_greeks;
    }

}

#endif //BLACKSCHOLES_GREEKS_H
