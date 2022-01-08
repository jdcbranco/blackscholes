#ifndef BSM_SOLVER_ANALYTICAL_INTERNALS_H
#define BSM_SOLVER_ANALYTICAL_INTERNALS_H

#include "solver.h"

namespace bsm {
    namespace internals {
        template<typename T>
        T calculate_european_forward(pricing<T> const& p) {
            return p.S * exp(-p.q*p.tau) -  p.K * exp(-p.r*p.tau);
        }

        template<typename T>
        T calculate_d1(pricing<T> const& p) {
            return (log(p.S / p.K) + (p.r - p.q + p.sigma * p.sigma / 2) * p.tau) / (p.sigma * sqrt(p.tau));
        }

        template<typename T>
        T calculate_d2(pricing<T> const& p) {
            return (log(p.S / p.K) + (p.r - p.q - p.sigma * p.sigma / 2) * p.tau) / (p.sigma * sqrt(p.tau));
        }

        template<typename T>
        T calculate_european_call(pricing<T> const& p) {
            auto d1 = calculate_d1(p);
            auto d2 = calculate_d2(p);
            return exp(-p.q * p.tau) * p.S * cdf<T>(d1) - exp(-p.r * p.tau) * p.K * cdf<T>(d2);
        }

        template<typename T>
        T calculate_european_put(pricing<T> const& p) {
            auto d1 = calculate_d1(p);
            auto d2 = calculate_d2(p);
            return exp(-p.r * p.tau) * p.K * cdf<T>(-d2) - exp(-p.q * p.tau) * p.S * cdf<T>(-d1);
        }

        template<typename T>
        T calculate_gamma(pricing<T> const& p) {
            auto const d1 = calculate_d1<T>(p);
            return exp(-p.q*p.tau)*pdf<T>(d1)/(p.S*p.sigma*sqrt(p.tau));
        }

        template<typename T>
        T calculate_vega(pricing<T> const& p) {
            auto const d1 = calculate_d1<T>(p);
            return p.S*exp(-p.q*p.tau)*pdf<T>(d1)*sqrt(p.tau);
        }
        template<typename T>
        T calculate_theta(pricing<T> const& p, double sign) {
            auto d1 = calculate_d1(p);
            auto d2 = calculate_d2(p);
            auto& S = p.S;
            auto& K = p.K;
            auto& sigma = p.sigma;
            auto& r = p.r;
            auto& q = p.q;
            auto& tau = p.tau;
            return - r * K * exp(-r*tau) * cdf<T>(d2*sign) * sign
            + q * S * exp(-q*tau) * cdf<T>(d1*sign) * sign
            - 0.5 * sigma * S * exp(-q * tau) * pdf<T>(d1) / sqrt(tau);
        }
    }
}


#endif //BSM_SOLVER_ANALYTICAL_INTERNALS_H
