#include <iostream>
#include "bsm.h"

var phi(var const& x) {
    return 0.5 *  (1+erf(x * one_div_root_two));
}

var european(CP cp, double K, var const& S, var const& variance, var const& tau, var const& r, var const& q) {
    var sigma = sqrt(variance);
    auto const d1 = (log(S / K) + (r - q + variance / 2) * tau) / (sigma * sqrt(tau));
    auto const d2 = (log(S / K) + (r - q - variance / 2) * tau) / (sigma * sqrt(tau));
    switch (cp) {
        case CP::call:
            return exp(-q * tau) * S * phi(d1) - exp(-r * tau) * K * phi(d2);
        case CP::put:
            return exp(-r * tau) * K * phi(-d2) - exp(-q * tau) * S * phi(-d1);
    }
}

var forward(double K, var const& S, var const& tau, var const& r, var const& q) {
    return S * exp(-q*tau) -  K * exp(- r*tau);
}

var forward(var const& S, var const& tau, var const& r, var const& q) {
    return S * exp ((r-q) * tau);
}

var implied_dividend(double K, var const& F, var const& S, var const& tau, var const& r, var q) {
    var f = abs(F - forward(K, S, tau, r, q));

    for(int i = 0; i<100; i++) {
        if(f<1e-9) {
            //std::cout << "Newton-Raphson for IDiv took " << i << " iterations" << std::endl;
            break;
        }
        auto [dfdq] = derivativesx(f, wrt(q));
        q.update(val(q-f/dfdq));
        f.update();
    }

    return q;
}

var implied_volatility(CP cp, double K, var const& P, var const& S, var const& tau, var const& r, var const& q, var sigma) {
    var variance = sigma*sigma;
    var f = abs(P - european(cp,K,S,variance,tau,r,q));

    for(int i = 0; i<100; i++) {
        if(f<1e-9) {
            //std::cout << "Newton-Raphson for IVol took " << i << " iterations" << std::endl;
            break;
        }
        auto [vega] = derivativesx(f, wrt(sigma));
        sigma.update(val(sigma-f/vega));
        variance.update();
        f.update();
    }

    return sigma;
}