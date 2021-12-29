#ifndef BLACKSCHOLES_ANALYTICAL_H
#define BLACKSCHOLES_ANALYTICAL_H

#include <array>
#include <autodiff/reverse/var.hpp>
#include <autodiff/reverse/var/eigen.hpp>
#include "options.h"
#include "greeks.h"

namespace analytical {
    using namespace autodiff;
    using namespace greeks;
    using namespace options;

    template<typename T>
    class analyticalmethod {
        market_variables<T> params_ptr;
    public:
        analyticalmethod(params<T> const& parameters): params_ptr{std::make_shared<params<T>>(parameters)} {}

        template<typename I> requires options::ClosedFormSolution<I,T>
        auto solve(I instrument) {
            return instrument.closed_form_price(*params_ptr);
        }
    };

    template<>
    class analyticalmethod<var> {
        market_variables<var> params_ptr;
    public:
        explicit analyticalmethod(params<var> const& parameters): params_ptr{std::make_shared<params<var>>(parameters)} {}

        template<typename I> requires ClosedFormSolution<I,var>
        auto solveForGreeks(I instrument) {
            var price = instrument.closed_form_price(*params_ptr);
            return std::make_tuple(price,greeks(price)); //<var,std::array<var,7>>
        }

        template<typename I> requires ClosedFormSolution<I,var>
        pricing<var> solve(I instrument) {
            auto price = instrument.closed_form_price(*params_ptr);
            pricing<var> pricing{price,params_ptr};
            return pricing;
        }

        std::array<var,6> greeks(var const& price) {
            auto& [S, sigma, tau, r, q] = *params_ptr;
            auto [delta, vega, theta, rho, psi] =  derivativesx(price, wrt(S, sigma, tau, r, q));
            auto [gamma] = derivativesx(delta, wrt(S));
            std::array<var,6> greeks{delta, gamma, vega, theta, rho, psi};
            return greeks;
        }

        template<typename I> requires ClosedFormSolution<I,var>
        auto imply_dividends(I instrument, var price, var q = 0.0) {
            params<var> local_params{params_ptr->S, params_ptr->sigma, params_ptr->tau, params_ptr->r, q};
            var f = abs(price - instrument.closed_form_price(local_params));

            for(int i = 0; i<100; i++) {
                if(f<1e-9) {
                    //std::cout << "Newton-Raphson for IDiv took " << i << " iterations" << std::endl;
                    break;
                }
                auto [dfdq] = derivativesx(f, wrt(q));
                q.update(val(q-f/dfdq));
                f.update();
            }

            return q;
        }

        template<typename I> requires ClosedFormSolution<I,var>
        auto imply_volatility(I instrument, var price, var sigma = 0.10) {
            params<var> local_params{params_ptr->S, sigma, params_ptr->tau, params_ptr->r, params_ptr->q};
            var f = abs(price - instrument.closed_form_price(local_params));

            for(int i = 0; i<100; i++) {
                if(f<1e-9) {
                    //std::cout << "Newton-Raphson for IVol took " << i << " iterations" << std::endl;
                    break;
                }
                auto [vega] = derivativesx(f, wrt(sigma));
                sigma.update(val(sigma-f/vega));
                f.update();
            }

            return sigma;
        }

    };

}


#endif //BLACKSCHOLES_ANALYTICAL_H
