#include "solver.h"

#include <optional>

#include <autodiff/reverse/var.hpp>
using namespace autodiff;

namespace bsm {

    struct var_pricing {
        var S;
        double K;
        var sigma;
        var tau;
        var r;
        var q;
        //Constructors
        var_pricing(forward const& instrument, mkt_params<double> mp):
        S{mp.S}, K{instrument.K}, sigma{mp.sigma},
        tau{static_cast<double>(time_between(mp.t,instrument.maturity).count())},
        r{mp.r}, q{mp.q} {}
        var_pricing(european const& instrument, mkt_params<double> mp):
        S{mp.S}, K{instrument.K}, sigma{mp.sigma},
        tau{static_cast<double>(time_between(mp.t,instrument.maturity).count())},
        r{mp.r}, q{mp.q} {}
        var_pricing(var_pricing const&) = default;
        var_pricing(var_pricing &&) noexcept = default;
    };

    var calculate_forward_var(var_pricing const& p) {
        return p.S * exp(-p.q*p.tau) -  p.K * exp(-p.r*p.tau);
    }

    var calculate_european_call_var(var_pricing const& p) {
        auto d1 = (log(p.S / p.K) + (p.r - p.q + p.sigma * p.sigma / 2) * p.tau) / (p.sigma * sqrt(p.tau));
        auto  d2 = (log(p.S / p.K) + (p.r - p.q - p.sigma * p.sigma / 2) * p.tau) / (p.sigma * sqrt(p.tau));
        return exp(-p.q * p.tau) * p.S * phi<var>(d1) - exp(-p.r * p.tau) * p.K * phi<var>(d2);
    }

    var calculate_european_put_var(var_pricing const& p) {
        auto d1 = (log(p.S / p.K) + (p.r - p.q + p.sigma*p.sigma / 2) * p.tau) / (p.sigma * sqrt(p.tau));
        auto d2 = (log(p.S / p.K) + (p.r - p.q - p.sigma*p.sigma / 2) * p.tau) / (p.sigma * sqrt(p.tau));
        return exp(-p.r * p.tau) * p.K * phi<var>(-d2) - exp(-p.q * p.tau) * p.S * phi<var>(-d1);
    }
    struct f_pricing_var: pricing {
        var_pricing vp;
        f_pricing_var(forward const& instrument, mkt_params<double> mp):
                pricing{instrument,mp},
                vp{instrument, mp}
        {}

        double price() override {
            return val(init_price());
        }

        double delta() override {
            return val(init_delta());
        }

        double gamma() override {
            return val(derivativesx(init_delta(), wrt(vp.S))[0]);
        }

        double vega() override {
            return val(derivativesx(init_price(), wrt(vp.sigma))[0]);
        }

        double theta() override {
            return val(derivativesx(init_price(), wrt(vp.tau))[0]);
        }

        double rho() override {
            return val(derivativesx(init_price(), wrt(vp.r))[0]);
        }

        double psi() override {
            return val(derivativesx(init_price(), wrt(vp.q))[0]);
        }
    protected:
        std::optional<var> price_;
        std::optional<var> delta_;
        var& init_price() {
            if(!price_) {
                price_ = calculate_forward_var(vp);
            }
            return price_.value();
        };

        var& init_delta() {
            if(!delta_) {
                delta_ = derivativesx(init_price(), wrt(vp.S))[0];
            }
            return delta_.value();
        }
    };

    struct e_pricing_var: pricing {
        var_pricing vp;
        e_pricing_var(european const& instrument, mkt_params<double> mp):
                pricing{instrument,mp},
                vp{instrument, mp}
        {}

        double price() override {
            return val(init_price());
        }

        double delta() override {
            return val(init_delta());
        }

        double gamma() override {
            return val(derivativesx(init_delta(), wrt(vp.S))[0]);
        }

        double vega() override {
            return val(derivativesx(init_price(), wrt(vp.sigma))[0]);
        }

        double theta() override {
            return val(derivativesx(init_price(), wrt(vp.tau))[0]);
        }

        double rho() override {
            return val(derivativesx(init_price(), wrt(vp.r))[0]);
        }

        double psi() override {
            return val(derivativesx(init_price(), wrt(vp.q))[0]);
        }
    protected:
        std::optional<var> price_;
        std::optional<var> delta_;
        virtual var& init_price() = 0;
        var& init_delta() {
            if(!delta_) {
                delta_ = derivativesx(init_price(), wrt(vp.S))[0];
            }
            return delta_.value();
        }
    };

    struct ec_pricing_var: e_pricing_var {
        ec_pricing_var(european_call const& instrument, mkt_params<double> mp): e_pricing_var{instrument, mp} {}
    protected:
        var& init_price() override {
            if(!price_) {
                price_ = calculate_european_call_var(vp);
            }
            return price_.value();
        }
    };

    struct ep_pricing_var: e_pricing_var {
        ep_pricing_var(european_put const& instrument, mkt_params<double> mp): e_pricing_var{instrument, mp} {}
    protected:
        var& init_price() override {
            if(!price_) {
                price_ = calculate_european_put_var(vp);
            }
            return price_.value();
        }
    };

    template<>
    std::unique_ptr<pricing> analytical_solver<autodiff_var>::operator()(forward& instrument) {
        f_pricing_var fp{instrument, mktParams};
        return std::make_unique<f_pricing_var>(fp);
    }

    template<>
    std::unique_ptr<pricing> analytical_solver<autodiff_var>::operator()(european_call& instrument) {
        ec_pricing_var fp{instrument, mktParams};
        return std::make_unique<ec_pricing_var>(fp);
    }

    template<>
    std::unique_ptr<pricing> analytical_solver<autodiff_var>::operator()(european_put& instrument) {
        ep_pricing_var fp{instrument, mktParams};
        return std::make_unique<ep_pricing_var>(fp);
    }

}
