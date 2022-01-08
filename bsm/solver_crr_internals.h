#ifndef BSM_SOLVER_CRR_INTERNALS_H
#define BSM_SOLVER_CRR_INTERNALS_H

#include "bintree.h"
#include "common.h"
#include "instruments.h"
#include "solver.h"

#include <vector>
#include <iostream>
#include <algorithm>
#include <execution>

namespace bsm {
    namespace internals {

        template<typename T>
        using calc_payoff_type = T(T const&);

        template<typename T>
        struct pricing_params {
            T S;
            double K;
            T sigma;
            T tau;
            T r;
            T q;
            pricing_params(instrument const& instrument, mkt_params<double> mp):
                    S{mp.S}, K{instrument.K}, sigma{mp.sigma},
                    tau{static_cast<double>(time_between(mp.t,instrument.maturity).count())},
                    r{mp.r}, q{mp.q} {}
            pricing_params(pricing_params const&) = default;
            pricing_params(pricing_params &&) noexcept = default;
        };

        template<typename T>
        struct generic_crr_pricing_method { //: pricing {
        protected:
            const int steps_;
            bintree<T> underlying_tree;
            bintree<T> premium_tree;
            T u_, d_, p_, discount_factor_;
        public:
            pricing_params<T> pp;
            generic_crr_pricing_method(instrument const& instrument, mkt_params<double> mp, int steps):
                    pp{instrument, mp},
                    underlying_tree{steps + 1}, premium_tree{steps + 1}, steps_{steps} {
                generate_underlying_tree();
            }
            generic_crr_pricing_method(instrument const& instrument, pricing_params<T> pp, int steps):
                    pp{pp},
                    underlying_tree{steps + 1}, premium_tree{steps + 1}, steps_{steps} {
                generate_underlying_tree();
            }

            void generate_underlying_tree() {
                int steps = premium_tree.size()-1;
                auto dt = pp.tau / steps;
                auto sqrt_dt = sqrt(dt);
                auto u = u_ = exp(pp.sigma * sqrt_dt);
                auto d = d_ = exp(-pp.sigma * sqrt_dt);
                p_ = (exp((pp.r - pp.q) * dt) - d) / (u - d);
                discount_factor_ = exp(-pp.r * dt);
                auto S = pp.S;
                underlying_tree.set(0, 0, S); //TODO Implement some syntatic sugar: underlying_tree(0,0) = S;
                std::vector<int> indices(steps+1);
                std::iota(indices.begin(),indices.end(), 0);
                for (int t = 1; t <= steps; ++t) {
                    auto start = indices.begin();
                    auto end = start+t+1;
                    auto [output, _ ] = underlying_tree(t);
                    transform(std::execution::par_unseq, start, end, output, [&S,&u,&d,&t](int i) {
                        return S*pow(u,t-i)*pow(d,i);
                    });
                }
            }

            T price() {
                return premium_tree.root();
            }

            T delta() {
                auto [V_start, V_end] = premium_tree(1);
                auto [S_start, S_end] = underlying_tree(1);
                auto V_up = *V_start;
                auto V_down = *(V_start+1);
                auto S_up = *S_start;
                auto S_down = *(S_start+1);
                return (V_up - V_down)/(S_up - S_down);
            }

            T gamma() {
                auto [V_start, V_end] = premium_tree(2);
                auto [S_start, S_end] = underlying_tree(2);
                auto V_uu = *V_start;
                auto V_ud = *(V_start+1);
                auto V_dd = *(V_start+2);
                auto S_uu = *S_start;
                auto S_ud = *(S_start+1);
                auto S_dd = *(S_start+2);
                return ( (V_uu-V_ud)/(S_uu-S_ud) - (V_ud-V_dd)/(S_ud-S_dd) )/((S_uu - S_dd)/2.0);
            }

            T theta() {
                auto V = premium_tree.root();
                auto V_ud = premium_tree(2,1);
                int steps = premium_tree.size()-1;
                auto dt = pp.tau / steps;
                return (V_ud-V)/(2*dt);
            }

            void solve(std::function<calc_payoff_type<T>> calc_payoff, bool early_exercise_possible) {
                auto steps = premium_tree.size();
                auto last_t = steps-1;
                auto p = p_;
                auto discount_factor = discount_factor_;
                {
                    auto [start, end] = underlying_tree(last_t);
                    auto [output, _] = premium_tree(last_t);
                    std::transform(std::execution::par_unseq, start, end, output, calc_payoff);
                }

                std::vector<int> indices(steps);
                std::iota(indices.begin(),indices.end(), 0);
                for(int t = last_t-1; t>=0; t--) {
                    auto start = indices.begin();
                    auto end = start+t+1;

                    auto [output, _] = premium_tree(t);
                    auto premium_next_step = get<0>(premium_tree(t+1));

                    std::transform(std::execution::par_unseq, start, end, output,
                              [premium_next_step, p, discount_factor, early_exercise_possible, &calc_payoff, t, this](int i) {
                                  auto premium_up = *(premium_next_step+i);
                                  auto premium_down = *(premium_next_step+i+1);
                                  auto continuation = (p*premium_up + (1.0-p)*premium_down)*discount_factor;
                                  return early_exercise_possible? std::max(continuation, calc_payoff(this->underlying_tree(t,i))) : continuation;
                              });
                }
            }

            auto underlying() const {
                return underlying_tree;
            }

            auto premium() const {
                return premium_tree;
            }

        };

        template<typename T>
        std::ostream& operator<<(std::ostream& out, generic_crr_pricing_method<T> const& crrtree) {
            out << crrtree.underlying();
            return out;
        }
    }
}


#endif //BSM_SOLVER_CRR_INTERNALS_H
