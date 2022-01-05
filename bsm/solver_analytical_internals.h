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
        T calculate_european_call(pricing<T> const& p) {
            auto d1 = (log(p.S / p.K) + (p.r - p.q + p.sigma * p.sigma / 2) * p.tau) / (p.sigma * sqrt(p.tau));
            auto d2 = (log(p.S / p.K) + (p.r - p.q - p.sigma * p.sigma / 2) * p.tau) / (p.sigma * sqrt(p.tau));
            return exp(-p.q * p.tau) * p.S * cdf<T>(d1) - exp(-p.r * p.tau) * p.K * cdf<T>(d2);
        }

        template<typename T>
        T calculate_european_put(pricing<T> const& p) {
            auto d1 = (log(p.S / p.K) + (p.r - p.q + p.sigma*p.sigma / 2) * p.tau) / (p.sigma * sqrt(p.tau));
            auto d2 = (log(p.S / p.K) + (p.r - p.q - p.sigma*p.sigma / 2) * p.tau) / (p.sigma * sqrt(p.tau));
            return exp(-p.r * p.tau) * p.K * cdf<T>(-d2) - exp(-p.q * p.tau) * p.S * cdf<T>(-d1);
        }

        template<typename T>
        T calculate_d1(pricing<T> const& p) {
            return (log(p.S / p.K) + (p.r - p.q + p.sigma*p.sigma / 2) * p.tau) / (p.sigma * sqrt(p.tau));
        }

        template<typename T>
        T calculate_d2(pricing<T> const& p) {
            return (log(p.S / p.K) + (p.r - p.q - p.sigma*p.sigma / 2) * p.tau) / (p.sigma * sqrt(p.tau));
        }
    }
}


#endif //BSM_SOLVER_ANALYTICAL_INTERNALS_H
