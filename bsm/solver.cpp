#include "solver.h"

#include <memory>
#include <cmath>

namespace bsm {

    double calculate_forward_d(pricing const& p) {
        return  p.S * exp(-p.q*p.tau) -  p.K * exp(-p.r*p.tau);
    }

    double calculate_european_call_d(pricing const& p) {
        auto const d1 = (log(p.S / p.K) + (p.r - p.q + p.sigma * p.sigma / 2) * p.tau) / (p.sigma * sqrt(p.tau));
        auto const d2 = (log(p.S / p.K) + (p.r - p.q - p.sigma * p.sigma / 2) * p.tau) / (p.sigma * sqrt(p.tau));
        return exp(-p.q * p.tau) * p.S * phi(d1) - exp(-p.r * p.tau) * p.K * phi(d2);
    }

    double calculate_european_put_d(pricing const& p) {
        auto const d1 = (log(p.S / p.K) + (p.r - p.q + p.sigma*p.sigma / 2) * p.tau) / (p.sigma * sqrt(p.tau));
        auto const d2 = (log(p.S / p.K) + (p.r - p.q - p.sigma*p.sigma / 2) * p.tau) / (p.sigma * sqrt(p.tau));
        return exp(-p.r * p.tau) * p.K * phi(-d2) - exp(-p.q * p.tau) * p.S * phi(-d1);
    }

    struct forward_pricing: pricing {
        forward_pricing(forward const& instrument, mkt_params<double> mp): pricing{instrument,mp} {}

        double price() override {
            return calculate_forward_d(*this);
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

    struct ec_pricing: pricing {
        ec_pricing(european_call const& instrument, mkt_params<double> mp): pricing{instrument,mp} {}

        double price() override {
            return calculate_european_call_d(*this);
        }

        double delta() override {
            auto const d1 = (log(S / K) + (r - q + sigma * sigma / 2) * tau) / (sigma * sqrt(tau));
            return exp(-q * tau)*phi(d1);
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

    struct ep_pricing: pricing {
        ep_pricing(european_put const& instrument, mkt_params<double> mp): pricing{instrument,mp} {}

        double price() override {
            return calculate_european_put_d(*this);
        }

        double delta() override {
            auto const d1 = (log(S / K) + (r - q + sigma * sigma / 2) * tau) / (sigma * sqrt(tau));
            return -exp(-q * tau)*phi(-d1);
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
    std::unique_ptr<pricing> analytical_solver<autodiff_off>::operator()(forward& instrument) {
        forward_pricing fp{instrument, mktParams};
        return std::make_unique<forward_pricing>(fp);
    }

    template<>
    std::unique_ptr<pricing> analytical_solver<autodiff_off>::operator()(european_call& instrument) {
        ec_pricing fp{instrument, mktParams};
        return std::make_unique<ec_pricing>(fp);
    }

    template<>
    std::unique_ptr<pricing> analytical_solver<autodiff_off>::operator()(european_put& instrument) {
        ep_pricing fp{instrument, mktParams};
        return std::make_unique<ep_pricing>(fp);
    }

}
