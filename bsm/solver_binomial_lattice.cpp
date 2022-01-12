
#include "solver.h"
#include "solver_lattice_internals.h"
#include "solver_american_internals.h"

#include <optional>

using namespace bsm::internals;

namespace bsm {

    struct sbl_method: american_method {

        pricing<long double> p;
        superpositioned_binomial_lattice_method<long double> sbl;

        sbl_method(american_put const& instrument, mkt_params<double> const& mp, int steps): sbl{instrument, steps}, p{instrument,mp} {
            sbl.solve(p);
        }
        sbl_method(sbl_method const&) = default;
        sbl_method(sbl_method &&) noexcept = default;

        double price() override {
            return sbl.price();
        }

        double delta() override {
            return sbl.delta();
        }

        double gamma() override {
            return sbl.gamma();
        }

        double vega() override {
            return sbl.vega();
        }

        double theta() override {
            return sbl.theta();
        }

        double rho() override {
            return sbl.rho();
        }

        double psi() override {
            return sbl.psi();
        }

        long double exercise_boundary(long double _tau) override {
            return sbl.exercise_boundary(_tau);
        }
    };

    template<>
    std::unique_ptr<american_method> sbl_solver<autodiff_off>::operator()(american_put& instrument) {
        sbl_method gp{instrument, mktParams, steps};
        return std::make_unique<sbl_method>(gp);
    }

}

