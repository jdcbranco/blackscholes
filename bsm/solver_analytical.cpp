#include "solver.h"
#include "solver_analytical_internals.h"

#include <memory>
#include <cmath>

using namespace bsm::internals;

namespace bsm {

    struct ef_analytical_pricing_method: pricing<double>, method {
        ef_analytical_pricing_method(european_forward const& instrument, mkt_params<double> mp): pricing{instrument, mp} {}

        double price() override {
            return calculate_european_forward<double>(*this);
        }

        double delta() override {
            return exp(-q*tau);
        }

        double gamma() override {
            return 0.0;
        }

        double vega() override {
            return 0.0;
        }

        double theta() override {
            return r*K*exp(-r*tau) - q*S*exp(-q*tau);
        }

        double rho() override {
            return tau*K*exp(-r*tau);
        }

        double psi() override {
            return -q*S*exp(-q*tau);
        }
    };

    struct ec_analytical_pricing_method: pricing<double>, method {
        ec_analytical_pricing_method(european_call const& instrument, mkt_params<double> mp): pricing{instrument, mp} {}

        double price() override {
            return calculate_european_call<double>(*this);
        }

        double delta() override {
            auto const d1 = (log(S / K) + (r - q + sigma * sigma / 2) * tau) / (sigma * sqrt(tau));
            return exp(-q * tau) * cdf(d1);
        }

        double gamma() override {
            auto const d1 = (log(S / K) + (r - q + sigma * sigma / 2) * tau) / (sigma * sqrt(tau));
            return exp(-q*tau)*exp(-0.5*d1*d1)/(S*sigma*sqrt(tau)*sqrt(2.0*std::numbers::pi));
        }

        //The following are not implemented yet
        double vega() override {
            return 0;
        }

        double theta() override {
            return 0;
        }

        double rho() override {
            return 0;
        }

        double psi() override {
            return 0;
        }
    };

    struct ep_analytical_pricing_method: pricing<double>, method {
        ep_analytical_pricing_method(european_put const& instrument, mkt_params<double> mp): pricing{instrument, mp} {}

        double price() override {
            return calculate_european_put<double>(*this);
        }

        double delta() override {
            auto const d1 = (log(S / K) + (r - q + sigma * sigma / 2) * tau) / (sigma * sqrt(tau));
            return -exp(-q * tau) * cdf(-d1);
        }

        double gamma() override {
            auto const d1 = (log(S / K) + (r - q + sigma * sigma / 2) * tau) / (sigma * sqrt(tau));
            return exp(-q*tau)*exp(-0.5*d1*d1)/(S*sigma*sqrt(tau)*sqrt(2.0*std::numbers::pi));
        }

        //The following are not implemented yet
        double vega() override {
            return 0;
        }

        double theta() override {
            return 0;
        }

        double rho() override {
            return 0;
        }

        double psi() override {
            return 0;
        }
    };

    template<>
    std::unique_ptr<method> analytical_solver<autodiff_off>::operator()(european_forward& instrument) {
        ef_analytical_pricing_method fp{instrument, mktParams};
        return std::make_unique<ef_analytical_pricing_method>(fp);
    }

    template<>
    std::unique_ptr<method> analytical_solver<autodiff_off>::operator()(european_call& instrument) {
        ec_analytical_pricing_method fp{instrument, mktParams};
        return std::make_unique<ec_analytical_pricing_method>(fp);
    }

    template<>
    std::unique_ptr<method> analytical_solver<autodiff_off>::operator()(european_put& instrument) {
        ep_analytical_pricing_method fp{instrument, mktParams};
        return std::make_unique<ep_analytical_pricing_method>(fp);
    }

}
