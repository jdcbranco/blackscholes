#include <iostream>
#include <autodiff/reverse/var.hpp>
#include <execution>
#include "crr.h"
#include "analytical.h"
#include "greeks.h"
#include "lsmc.h"

using namespace autodiff;

int main() {
    using namespace greeks;

    using std::cout;
    using std::endl;

    //Example of European Call pricing and greek calculation using autodiff
    double const K = 15000.0;

    options::params<var> main_params{15000.0, 0.21, 0.5, 0.0, 0.0};
    auto&& [S, sigma, tau, r, q] = main_params;

    options::european_call<var> europeanCall{K};
    analytical::analyticalmethod<var> bsm_reverse_mode{main_params};

    auto call = bsm_reverse_mode.solve(europeanCall);
    auto [delta, gamma, vega, theta, rho, psi] = all_greeks(call); // greeks;

    std::cout << "call = " << call << std::endl;
    std::cout << "delta = " << delta << std::endl;
    std::cout << "gamma = " << gamma << std::endl;
    std::cout << "vega = " << vega << std::endl;
    std::cout << "theta = " << theta << std::endl;
    std::cout << "rho = " << rho << std::endl;
    std::cout << "psi = " << psi << std::endl;
    std::cout << "\n\n";

    //Example of Implied Dividend calculation using Newton-Raphson method and autodiff
    options::forward<var> fwdInstrument{K};
    auto fwd = bsm_reverse_mode.solve(fwdInstrument);
    auto impl_div = bsm_reverse_mode.imply_dividends(fwdInstrument, (var)fwd);

    std::cout << "fwd = " << fwd << std::endl;
    std::cout << std::boolalpha;
    std::cout << "impl div = " << impl_div << " which is " << (abs(impl_div - q) <= 1.0e-4) << std::endl;

    //Example of Implied Volatility calculation using Newton-Raphson method and autodiff
    auto iv = bsm_reverse_mode.imply_volatility(europeanCall, val((var)call));

    std::cout << "impl vol = " << iv << " which is " << (abs(iv - sigma) <= 1.0e-9) << std::endl;

    auto K_ = K;
    auto S_ = val(S);
    auto sigma_ = val(sigma);
    auto tau_ = val(tau);
    auto r_ = val(r);
    auto q_ = val(q);

    options::params<double> american_params{S_, sigma_, tau_, r_, q_};

    {
        bool run_lsmc = false;

        if (run_lsmc) {
            cout << "Least Squares Monte Carlo using doubles:" << endl;
            options::american_call<double> americanCall{K_};
            lsmc::lsmc_method<double> lsmcMethod{american_params};

            int total_sims = 100;
            std::vector<int> mc_input(total_sims);
            std::vector<double> mc_simulations(total_sims);

            std::transform(std::execution::par, mc_input.begin(), mc_input.end(), mc_simulations.begin(),
                           [&lsmcMethod, &americanCall](int i) { //K_,S_,tau_,r_,q_
                               options::pricing<double> mc_pricing = lsmcMethod.solve(americanCall, 180, 1000);
                               return mc_pricing.price();
                               //return american(CP::call, K_, S_, variance_, tau_, r_, q_, 180, 1000);
                           });
            auto american_call_price =
                    std::reduce(std::execution::par, mc_simulations.begin(), mc_simulations.end(), 0.0) / total_sims;

            std::cout << "american call (lsmc) = " << american_call_price << std::endl;
        }
    }

    {
        //This section is too slow to be worthwhile
        bool run_lsmc = false;

        if (run_lsmc) {
            cout << "Least Squares Monte Carlo using Autodiff:" << endl;
            options::params<var> american_params{S_, sigma_, tau_, r_, q_};
            options::american_call<var> americanCall{K_};
            lsmc::lsmc_method<var> lsmcMethod{american_params};

            int total_sims = 1;
            std::vector<int> mc_input(total_sims);
            std::vector<double> mc_simulations(total_sims);

            std::transform(std::execution::par, mc_input.begin(), mc_input.end(), mc_simulations.begin(),
                           [&lsmcMethod, &americanCall](int i) { //K_,S_,tau_,r_,q_
                               options::pricing<var> mc_pricing = lsmcMethod.solve(americanCall, 180, 10);
                               return val(mc_pricing.price());
                               //return american(CP::call, K_, S_, variance_, tau_, r_, q_, 180, 1000);
                           });
            auto american_call_price =
                    std::reduce(std::execution::par, mc_simulations.begin(), mc_simulations.end(), 0.0 ) / total_sims;

            std::cout << "american call (lsmc) = " << american_call_price << std::endl;
        }
    }

    {
        auto binomial_steps = 400;
        cout << "Binomial Tree Example using doubles: " << endl;

        crr::crrmethod<double> crr{american_params, binomial_steps};
        options::european_call<double> europeanCall2{K_};
        auto result = crr.solve(europeanCall2);

        if (binomial_steps <= 5) {
            cout << std::fixed << crr.underlying() << endl;
            cout << crr.premium() << endl;
        }

        cout << "european call (crr) = " << result << " with " << binomial_steps << " steps " << endl;
    }

    {
        //Too slow as well, uses too much memory
        auto binomial_steps = 10;
        cout << "Binomial Tree Example using Autodiff: " << endl;
        options::params<var> american_params{S_, sigma_, tau_, r_, q_};
        crr::crrmethod<var> crr{american_params, binomial_steps};
        options::european_call<var> europeanCall2{K_};
        auto result = crr.solve(europeanCall2);

        if (binomial_steps <= 5) {
            cout << std::fixed << crr.underlying() << endl;
            cout << crr.premium() << endl;
        }

        auto delta = greeks::delta(result);
        auto gamma = greeks::gamma(result);
        auto vega = greeks::vega(result);

        cout << "european call (crr) = " << result << " with " << binomial_steps << " steps " << endl;
        std::cout << "delta (crr) = " << delta << std::endl;
        std::cout << "gamma (crr) = " << gamma << std::endl;
        std::cout << "vega (crr) = " << vega << std::endl;
    }


    return 0;
}