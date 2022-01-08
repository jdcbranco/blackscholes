#include "solver.h"
#include "solver_analytical_internals.h"

#include <concepts>
#include <optional>
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
using namespace autodiff;
using namespace bsm::internals;

namespace bsm {

    struct dual_pricing_method: pricing<double>, method {

        dual_pricing_method(instrument const& inst, mkt_params<double> const& mp, std::function<pricing_function<dual2nd>> const& calc):
        pricing{inst,mp} , dp{inst,mp}, calc{calc} {
            auto [_price, _delta, _gamma] = derivatives(calc, wrt(dp.S, dp.S), at(dp));
            price_ = _price;
            delta_ = _delta;
            gamma_ = _gamma;
            vega_ = derivative(calc, wrt(dp.sigma), at(dp));
            theta_ = derivative(calc, wrt(dp.tau), at(dp));
            rho_ = derivative(calc, wrt(dp.r), at(dp));
            psi_ = derivative(calc, wrt(dp.q), at(dp));
        }

        double price() override {
            return val(price_);
        }

        double delta() override {
            return val(delta_);
        }

        double gamma() override {
            return val(gamma_);
        }

        double vega() override {
            return val(vega_);
        }

        double theta() override {
            return -val(theta_);
        }

        double rho() override {
            return val(rho_);
        }

        double psi() override {
            return val(psi_);
        }

    protected:
        pricing<dual2nd> dp;
        std::function<pricing_function<dual2nd>> calc;
        dual2nd price_;
        dual2nd delta_;
        dual2nd gamma_;
        dual2nd vega_;
        dual2nd theta_;
        dual2nd rho_;
        dual2nd psi_;
    };

    template<>
    std::unique_ptr<method> analytical_solver<autodiff_dual>::operator()(european_forward& instrument) {
        dual_pricing_method gp{instrument, mktParams, [](pricing<dual2nd> p) { return calculate_european_forward<dual2nd>(p); }};
        return std::make_unique<dual_pricing_method>(gp);
    }

    template<>
    std::unique_ptr<method> analytical_solver<autodiff_dual>::operator()(european_call& instrument) {
        dual_pricing_method gp{instrument, mktParams, [](pricing<dual2nd> p) { return calculate_european_call<dual2nd>(p); }};
        return std::make_unique<dual_pricing_method>(gp);
    }

    template<>
    std::unique_ptr<method> analytical_solver<autodiff_dual>::operator()(european_put& instrument) {
        dual_pricing_method gp{instrument, mktParams, [](pricing<dual2nd> p) { return calculate_european_put<dual2nd>(p); }};
        return std::make_unique<dual_pricing_method>(gp);
    }

}
