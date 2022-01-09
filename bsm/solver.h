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
    const double one_div_root_two_pi = 1.0/sqrt(2.0*std::numbers::pi);

    template<typename T = double>
    inline T cdf(T const& x)
    {
        return 0.5 * (1.0 + erf(x * one_div_root_two));
    }

    template<typename T = double>
    inline T pdf(T const& x) {
        return exp(-0.5*x*x)*one_div_root_two_pi;
    }

    struct method {
        virtual double price() = 0;
        virtual double delta() = 0;
        virtual double gamma() = 0;
        virtual double vega() = 0;
        virtual double theta() = 0;
        virtual double rho() = 0;
        virtual double psi() = 0;
    };

    struct american_method: method {
        virtual long double exercise_boundary(long double _tau) = 0;
    };

    template<typename T = double>
    struct pricing {
        T S;
        T K;
        T sigma;
        T tau;
        T r;
        T q;
        pricing(instrument const& instrument, mkt_params<double> const& mp):
                S{mp.S}, K{instrument.K}, sigma{mp.sigma},
                tau{static_cast<double>(time_between(mp.t,instrument.maturity).count())},
                r{mp.r}, q{mp.q} {}
        pricing(T const& S, T const& K, T const& sigma, T const& tau, T const& r, T const& q):
                S{S}, K{K}, sigma{sigma},
                tau{tau},
                r{r}, q{q} {}
        pricing(pricing const&) = default;
        pricing(pricing &&) noexcept = default;
        std::unique_ptr<pricing<T>> clone(T const& S, T const& tau) {
            pricing<T> copy{*this};
            copy.S = S;
            copy.tau = tau;
            return std::make_unique<pricing<T>>(copy);
        }
    };

    template<typename T>
    using pricing_function = T(pricing<T> const&);

//    template<typename T>
//    using exercise_boundary_function = T(T const&, T const&);

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

        std::unique_ptr<method> operator()(european_forward& instrument);
        std::unique_ptr<method> operator()(european_call& instrument);
        std::unique_ptr<method> operator()(european_put& instrument);
    };

    template<typename AD = autodiff_off>
    struct crr_solver {
        mkt_params<double> mktParams;
        const int steps;
        const int extra_steps;
    public:
        inline crr_solver(mkt_params<double> const& mktParams, int steps, int extra_steps = 0):
            mktParams{mktParams}, steps{steps}, extra_steps{extra_steps} {
            assert(("Extra steps must be even",extra_steps%2==0));
        }
        inline crr_solver(mkt_params<long double> const& mktParams, int steps, int extra_steps = 0): mktParams{mktParams}, steps{steps}, extra_steps{extra_steps} {}
        inline crr_solver(crr_solver const&) = default;
        inline crr_solver(crr_solver &&) noexcept = default;

        std::unique_ptr<method> operator()(european_forward& instrument);
        std::unique_ptr<method> operator()(european_call& instrument);
        std::unique_ptr<method> operator()(european_put& instrument);
        std::unique_ptr<method> operator()(american_call& instrument);
        std::unique_ptr<method> operator()(american_put& instrument);
    };

    template<typename AD = autodiff_off>
    struct qdplus_solver {
        mkt_params<long double> mktParams;
        inline qdplus_solver(mkt_params<long double> const& mktParams): mktParams{mktParams} {}
        inline qdplus_solver(qdplus_solver const&) = default;
        inline qdplus_solver(qdplus_solver &&) noexcept = default;

        std::unique_ptr<american_method> operator()(american_put& instrument);
        std::unique_ptr<american_method> operator()(american_call& instrument);

    };

    template<typename AD = autodiff_off>
    struct fastamerican_solver {
        mkt_params<double> mktParams;
        //method parameters
        const int l, m, n;
    public:
        inline fastamerican_solver(mkt_params<double> const& mktParams, int l, int m, int n): mktParams{mktParams}, l{l}, m{m}, n{n} {}
        inline fastamerican_solver(fastamerican_solver const&) = default;
        inline fastamerican_solver(fastamerican_solver &&) noexcept = default;

    };

}



#endif //BSM_SOLVER_H
