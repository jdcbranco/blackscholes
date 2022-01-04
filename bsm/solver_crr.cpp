#include "solver.h"

#include "solver_crr_internals.h"

using namespace bsm::internals;

namespace bsm {

    struct crr_pricing_method: pricing<double>, method {
        std::function<calc_payoff_type<double>> calc_payoff;
        generic_crr_pricing_method<double> crr;
        const instrument instrument_;
        const int steps;
        const bool early_exercise;

        crr_pricing_method(european const& instrument, mkt_params<double> mp, int steps):
        pricing{instrument,mp}, crr{instrument, mp, steps}, calc_payoff{[&instrument](double price) { return instrument.payoff(price); }}, steps{steps}, instrument_{instrument}, early_exercise{false}
        {
            crr.solve(calc_payoff, early_exercise);
        }

        crr_pricing_method(american const& instrument, mkt_params<double> mp, int steps):
                pricing{instrument,mp}, crr{instrument, mp, steps}, calc_payoff{[&instrument](double price) { return instrument.payoff(price); }}, steps{steps}, instrument_{instrument}, early_exercise{true}
        {
            crr.solve(calc_payoff, early_exercise);
        }

        double price() override {
            return crr.price();
        }

        double delta() override {
            return crr.delta();
        }

        double gamma() override {
            return crr.gamma();
        }

        double vega() override {
            pricing_params<double> bumped_up {crr.pp };
            bumped_up.sigma *= exp(0.01);
            generic_crr_pricing_method<double> bumped_up_crr{instrument_, bumped_up, steps};
            bumped_up_crr.solve(calc_payoff, early_exercise);
            return (bumped_up_crr.price() - crr.price()) / (bumped_up.sigma - crr.pp.sigma);
        }

        double theta() override {
            return crr.theta();
        }

        double rho() override {
            pricing_params<double> bumped_up {crr.pp };
            bumped_up.r *= exp(0.01);
            generic_crr_pricing_method<double> bumped_up_crr{instrument_, bumped_up, steps};
            bumped_up_crr.solve(calc_payoff,early_exercise);
            return (bumped_up_crr.price() - crr.price()) / (bumped_up.r - crr.pp.r);
        }

        double psi() override {
            pricing_params<double> bumped_up {crr.pp};
            bumped_up.q *= exp(0.01);
            generic_crr_pricing_method<double> bumped_up_crr{instrument_, bumped_up, steps};
            bumped_up_crr.solve(calc_payoff,early_exercise);
            return (bumped_up_crr.price() - crr.price()) / (bumped_up.q - crr.pp.q);
        }

    };

    template<>
    std::unique_ptr<method> crr_solver<autodiff_off>::operator()(european_forward& instrument) {
        crr_pricing_method gp{instrument, mktParams, steps};
        return std::make_unique<crr_pricing_method>(gp);
    }

    template<>
    std::unique_ptr<method> crr_solver<autodiff_off>::operator()(european_call& instrument) {
        crr_pricing_method gp{instrument, mktParams, steps };
        return std::make_unique<crr_pricing_method>(gp);
    }

    template<>
    std::unique_ptr<method> crr_solver<autodiff_off>::operator()(european_put& instrument) {
        crr_pricing_method gp{instrument, mktParams, steps };
        return std::make_unique<crr_pricing_method>(gp);
    }

    template<>
    std::unique_ptr<method> crr_solver<autodiff_off>::operator()(american_call& instrument) {
        crr_pricing_method gp{instrument, mktParams, steps };
        return std::make_unique<crr_pricing_method>(gp);
    }

    template<>
    std::unique_ptr<method> crr_solver<autodiff_off>::operator()(american_put& instrument) {
        crr_pricing_method gp{instrument, mktParams, steps };
        return std::make_unique<crr_pricing_method>(gp);
    }

}

