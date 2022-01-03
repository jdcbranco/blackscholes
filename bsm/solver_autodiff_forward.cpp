#include "solver.h"

#include <concepts>
#include <optional>

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
using namespace autodiff;

namespace bsm {

    template<typename T = dual>
    struct dual_pricing {
        T S;
        double K;
        T sigma;
        T tau;
        T r;
        T q;
        dual_pricing(instrument const& instrument, mkt_params<double> mp):
                S{mp.S}, K{instrument.K}, sigma{mp.sigma},
                tau{static_cast<double>(time_between(mp.t,instrument.maturity).count())},
                r{mp.r}, q{mp.q} {}
        dual_pricing(dual_pricing const&) = default;
        dual_pricing(dual_pricing &&) noexcept = default;
    };

    template<typename T>
    using calc_type = T(dual_pricing<T> const&);

    template<typename T = dual>
    T calculate_forward_dual(dual_pricing<T> const& p) {
        return p.S * exp(-p.q*p.tau) -  p.K * exp(-p.r*p.tau);
    }

    template<typename T = dual>
    T calculate_european_call_dual(dual_pricing<T> const& p) {
        auto d1 = (log(p.S / p.K) + (p.r - p.q + p.sigma * p.sigma / 2) * p.tau) / (p.sigma * sqrt(p.tau));
        auto  d2 = (log(p.S / p.K) + (p.r - p.q - p.sigma * p.sigma / 2) * p.tau) / (p.sigma * sqrt(p.tau));
        return exp(-p.q * p.tau) * p.S * phi<T>(d1) - exp(-p.r * p.tau) * p.K * phi<T>(d2);
    }

    template<typename T = dual>
    T calculate_european_put_dual(dual_pricing<T> const& p) {
        auto d1 = (log(p.S / p.K) + (p.r - p.q + p.sigma*p.sigma / 2) * p.tau) / (p.sigma * sqrt(p.tau));
        auto d2 = (log(p.S / p.K) + (p.r - p.q - p.sigma*p.sigma / 2) * p.tau) / (p.sigma * sqrt(p.tau));
        return exp(-p.r * p.tau) * p.K * phi<T>(-d2) - exp(-p.q * p.tau) * p.S * phi<T>(-d1);
    }

    struct dual_pricing_method: pricing {

        dual_pricing_method(instrument const& inst, mkt_params<double> const& mp, std::function<calc_type<dual2nd>> const& calc):
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
            return val(theta_);
        }

        double rho() override {
            return val(rho_);
        }

        double psi() override {
            return val(psi_);
        }

    protected:
        dual_pricing<dual2nd> dp;
        std::function<calc_type<dual2nd>> calc;
        dual2nd price_;
        dual2nd delta_;
        dual2nd gamma_;
        dual2nd vega_;
        dual2nd theta_;
        dual2nd rho_;
        dual2nd psi_;
    };

    template<>
    std::unique_ptr<pricing> analytical_solver<autodiff_dual>::operator()(forward& instrument) {
        dual_pricing_method gp{instrument, mktParams, [](dual_pricing<dual2nd> p) { return calculate_forward_dual(p); }};
        return std::make_unique<dual_pricing_method>(gp);
    }

    template<>
    std::unique_ptr<pricing> analytical_solver<autodiff_dual>::operator()(european_call& instrument) {
        dual_pricing_method gp{instrument, mktParams, [](dual_pricing<dual2nd> p) { return calculate_european_call_dual(p); }};
        return std::make_unique<dual_pricing_method>(gp);
    }

    template<>
    std::unique_ptr<pricing> analytical_solver<autodiff_dual>::operator()(european_put& instrument) {
        dual_pricing_method gp{instrument, mktParams, [](dual_pricing<dual2nd> p) { return calculate_european_put_dual(p); }};
        return std::make_unique<dual_pricing_method>(gp);
    }

}
