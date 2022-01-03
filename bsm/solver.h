#ifndef BSM_SOLVER_H
#define BSM_SOLVER_H

#include "common.h"
#include "instruments.h"
#include "bintree.h"

#include <concepts>
#include <iostream>
#include <sstream>
#include <memory>
#include <cmath>

namespace bsm {

    const double one_div_root_two = 1.0/sqrt(2.0);

    template<typename T = double>
    inline T phi(T const& x)
    {
        return 0.5 * (1.0 + erf(x * one_div_root_two));
    }

    struct pricing {
        double S;
        double K;
        double sigma;
        double tau;
        double r;
        double q;
        //Constructors
        pricing(instrument const& instrument, mkt_params<double> mp):
                S{mp.S}, K{instrument.K}, sigma{mp.sigma},
                tau{static_cast<double>(time_between(mp.t,instrument.maturity).count())},
                r{mp.r}, q{mp.q} {}
        pricing(forward const& instrument, mkt_params<double> mp):
            S{mp.S}, K{instrument.K}, sigma{mp.sigma},
            tau{static_cast<double>(time_between(mp.t,instrument.maturity).count())},
            r{mp.r}, q{mp.q} {}
        pricing(european const& instrument, mkt_params<double> mp):
            S{mp.S}, K{instrument.K}, sigma{mp.sigma},
            tau{static_cast<double>(time_between(mp.t,instrument.maturity).count())},
            r{mp.r}, q{mp.q} {}
        pricing(pricing const&) = default;
        pricing(pricing &&) noexcept = default;
        //Public methods
        virtual double price() = 0;
        virtual double delta() = 0;
        virtual double gamma() = 0;
        virtual double vega() = 0;
        virtual double theta() = 0;
        virtual double rho() = 0;
        virtual double psi() = 0;
    };

    class autodiff_off;
    class autodiff_dual;
    class autodiff_var;

    template<typename AD = autodiff_off>
    class analytical_solver {
        mkt_params<double> mktParams;
    public:
        inline analytical_solver(mkt_params<double> const& mktParams): mktParams{mktParams} {}
        inline analytical_solver(analytical_solver const&) = default;
        inline analytical_solver(analytical_solver &&) noexcept = default;

        std::unique_ptr<pricing> operator()(forward& instrument);
        std::unique_ptr<pricing> operator()(european_call& instrument);
        std::unique_ptr<pricing> operator()(european_put& instrument);
    };

    template<typename AD = autodiff_off>
    struct crr_solver {
        mkt_params<double> mktParams;
        const int steps;
    public:
        inline crr_solver(mkt_params<double> const& mktParams, int steps): mktParams{mktParams}, steps{steps} {}
        inline crr_solver(crr_solver const&) = default;
        inline crr_solver(crr_solver &&) noexcept = default;

        std::unique_ptr<pricing> operator()(forward& instrument);
        std::unique_ptr<pricing> operator()(european_call& instrument);
        std::unique_ptr<pricing> operator()(european_put& instrument);
        std::unique_ptr<pricing> operator()(american_call& instrument);
        std::unique_ptr<pricing> operator()(american_put& instrument);
    };

}



#endif //BSM_SOLVER_H
