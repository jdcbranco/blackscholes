#include "solver.h"
#include "solver_analytical_internals.h"

#include <optional>
#include <autodiff/reverse/var.hpp>

using namespace autodiff;
using namespace bsm::internals;

namespace bsm {

    using var_pricing = pricing<var>;

    struct var_pricing_method: pricing<double>, method {
        var_pricing vp;
        var_pricing_method(european const& instrument, mkt_params<double> mp, std::function<pricing_function<var>> const& calc):
                pricing{instrument,mp},
                vp{instrument, mp},
                calc{calc}
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
            return -val(derivativesx(init_price(), wrt(vp.tau))[0]);
        }

        double rho() override {
            return val(derivativesx(init_price(), wrt(vp.r))[0]);
        }

        double psi() override {
            return val(derivativesx(init_price(), wrt(vp.q))[0]);
        }
    protected:
        std::function<pricing_function<var>> calc;
        std::optional<var> price_;
        std::optional<var> delta_;
        var& init_price() {
            if(!price_) {
                price_ = calc(vp);
            }
            return price_.value();
        }
        var& init_delta() {
            if(!delta_) {
                delta_ = derivativesx(init_price(), wrt(vp.S))[0];
            }
            return delta_.value();
        }
    };

    template<>
    std::unique_ptr<method> analytical_solver<autodiff_var>::operator()(european_forward& instrument) {
        var_pricing_method fp{instrument, mktParams, calculate_european_forward<var>};
        return std::make_unique<var_pricing_method>(fp);
    }

    template<>
    std::unique_ptr<method> analytical_solver<autodiff_var>::operator()(european_call& instrument) {
        var_pricing_method fp{instrument, mktParams, calculate_european_call<var>};
        return std::make_unique<var_pricing_method>(fp);
    }

    template<>
    std::unique_ptr<method> analytical_solver<autodiff_var>::operator()(european_put& instrument) {
        var_pricing_method fp{instrument, mktParams, calculate_european_put<var>};
        return std::make_unique<var_pricing_method>(fp);
    }

}
