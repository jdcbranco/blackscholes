#include "solver.h"

#include "solver_crr_internals.h"
#include "solver_american_internals.h"

#include <optional>

using namespace bsm::internals;

namespace bsm {

    struct crr_pricing_method: pricing<double>, american_method {
        std::function<calc_payoff_type<double>> calc_payoff;
        generic_crr_pricing_method<double> crr;
        const instrument instrument_;
        const int steps;
        const bool early_exercise;
        std::optional<std::vector<double>> boundary;

        crr_pricing_method(european const& instrument, mkt_params<double> mp, int steps):
        pricing{instrument,mp}, crr{instrument, mp, steps}, calc_payoff{[&instrument](double price) { return instrument.payoff(price); }}, steps{steps}, instrument_{instrument}, early_exercise{false}
        {
            crr.solve(calc_payoff, early_exercise);
        }

        crr_pricing_method(american const& instrument, mkt_params<double> mp, int steps, int extra = 0):
                pricing{instrument,mp}, crr{instrument, mp, steps+extra, extra}, calc_payoff{[&instrument](double price) { return instrument.payoff(price); }}, steps{steps}, instrument_{instrument}, early_exercise{true}
        {
            boundary = crr.solve(calc_payoff, early_exercise);
            if(extra>0) {
                boundary->erase(boundary->begin(), boundary->begin()+extra);
            }
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
            if(bumped_up.q!=0)
                bumped_up.q *= exp(0.01);
            else
                bumped_up.q += 0.01;
            generic_crr_pricing_method<double> bumped_up_crr{instrument_, bumped_up, steps};
            bumped_up_crr.solve(calc_payoff,early_exercise);
            return (bumped_up_crr.price() - crr.price()) / (bumped_up.q - crr.pp.q);
        }

        long double exercise_boundary(long double _tau) override {
            bool call = instrument_.type==instrument_type::call;
            bool put = instrument_.type==instrument_type::put;
            if(never_optimal_exercise<double>(*this,instrument_.type)) {
                return call? INFINITY : 0.0;
            }
            if(tau==0) {
                return exercise_boundary_at_maturity<double>(*this,instrument_.type);
            }

            if(boundary) {
                std::cout << "Boundary size: "<< boundary.value().size() << std::endl;
                int index = steps*(1.0- _tau/tau);
                if(index >=0 && index <=steps) {
                    return boundary.value()[index];
                }
            }

            return NAN;
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
    std::unique_ptr<american_method> crr_solver<autodiff_off>::operator()(american_call& instrument) {
        crr_pricing_method gp{instrument, mktParams, steps, extra_steps };
        return std::make_unique<crr_pricing_method>(gp);
    }

    template<>
    std::unique_ptr<american_method> crr_solver<autodiff_off>::operator()(american_put& instrument) {
    crr_pricing_method gp{instrument, mktParams, steps, extra_steps };
        return std::make_unique<crr_pricing_method>(gp);
    }

}

