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
#include <utility>

namespace bsm {
    namespace internals {

        template<typename T>
        using calc_payoff_type = T(T const&);

        template<typename T>
        struct pricing_params {
            T S;
            long double K;
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
        struct generic_crr_pricing_method {
        protected:
            const int steps;
            const int shift;
            const instrument_type type;
            bintree<T> underlying_tree;
            bintree<std::pair<T,bool>> premium_tree;
            T u_, d_, p_, discount_factor_;
        public:
            pricing_params<T> pp;
            generic_crr_pricing_method(instrument const& instrument, mkt_params<double> mp, int steps, int shift = 0):
                    pp{instrument, mp},
                    underlying_tree{steps + 1}, premium_tree{steps + 1}, steps{steps}, shift{shift}, type{instrument.type} {
                generate_underlying_tree();
            }
            generic_crr_pricing_method(instrument const& instrument, pricing_params<T> pp, int steps, int shift = 0):
                    pp{pp},
                    underlying_tree{steps + 1}, premium_tree{steps + 1}, steps{steps}, shift{shift}, type{instrument.type} {
                generate_underlying_tree();
            }

            T pt(int i, int j) {
                return premium_tree(i + shift, j + shift / 2).first;
            }

            T ut(int i, int j) {
                return underlying_tree(i + shift, j + shift / 2);
            }

            void generate_underlying_tree() {
                auto dt = pp.tau / (steps - shift);
                auto sqrt_dt = sqrt(dt);
                auto u = u_ = exp(pp.sigma * sqrt_dt);
                auto d = d_ = exp(-pp.sigma * sqrt_dt);
                p_ = (exp((pp.r - pp.q) * dt) - d) / (u - d);
                discount_factor_ = exp(-pp.r * dt);
                auto S = pp.S;
                underlying_tree.set(0, 0, S);
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
                return pt(0,0);
            }

            T delta() {
                auto V_up   = pt(1,0);
                auto V_down = pt(1,1);
                auto S_up   = ut(1,0);
                auto S_down = ut(1,1);
                return (V_up - V_down)/(S_up - S_down);
            }

            T gamma() {
                auto V_uu = pt(2,0);
                auto V_ud = pt(2,1);
                auto V_dd = pt(2,2);
                auto S_uu = ut(2,0);
                auto S_ud = ut(2,1);
                auto S_dd = ut(2,2);
                return ( (V_uu-V_ud)/(S_uu-S_ud) - (V_ud-V_dd)/(S_ud-S_dd) )/((S_uu - S_dd)/2.0);
            }

            T theta() {
                auto V    = pt(0,0);
                auto V_ud = pt(2,1);
                auto dt = pp.tau / (steps - shift);
                return (V_ud-V)/(2*dt);
            }

            std::vector<T> solve(std::function<calc_payoff_type<T>> calc_payoff, bool early_exercise_possible) {
                auto last_t = steps;
                auto p = p_;
                auto discount_factor = discount_factor_;
                {
                    auto [start, end] = underlying_tree(last_t);
                    auto [output, _] = premium_tree(last_t);
                    std::transform(std::execution::par_unseq, start, end, output, [&calc_payoff](T const& price) {
                        return std::make_pair(calc_payoff(price),false);
                    }); //calc_payoff
                }

                std::vector<int> indices(premium_tree.size());
                std::vector<T> boundary(premium_tree.size());
                std::iota(indices.begin(),indices.end(), 0);
                for(int t = last_t-1; t>=0; t--) {
                    auto start = indices.begin();
                    auto end = start+t+1;

                    auto [output, _] = premium_tree(t);
                    auto premium_next_step = get<0>(premium_tree(t+1));

                    std::transform(std::execution::par_unseq, start, end, output,
                              [premium_next_step, p, discount_factor, early_exercise_possible, &calc_payoff, t, this](int i) {
                                  std::pair<T,bool> premium_up = *(premium_next_step+i);
                                  //auto premium_up = *(premium_next_step+i);
                                  std::pair<T,bool> premium_down = *(premium_next_step+i+1);
                                  auto continuation = (p*premium_up.first + (1.0-p)*premium_down.first)*discount_factor;
                                  if(early_exercise_possible) {
                                      T payoff = calc_payoff(this->underlying_tree(t,i));
                                      if(payoff > continuation) {
                                          return std::make_pair(payoff, true);
                                      } else {
                                          return std::make_pair(continuation,false);
                                      }
                                  } else {
                                      return std::make_pair(continuation,false);
                                  }
                                  //return early_exercise_possible? std::make_pair(std::max(continuation, calc_payoff(this->underlying_tree(t,i))),true) : std::make_pair(continuation,false);
                              });

                    if(early_exercise_possible) {
                        //calc exercise boundary
                        auto[start, end] = premium_tree(t);

                        int b = -1;
                        if (type == instrument_type::put) {
                            for (auto it = start; it != end; it++) {
                                if (it->second) {
                                    b = it - start;
                                    break;
                                }
                            }
                            if(b>0) {
                                //This approximation is based on paper "Discrete and continuous time approximations of the optiomal exercise boundary of American options - Basso, Nardon, Pianca"
                                auto den = premium_tree(t, b - 1).first - premium_tree(t, b).first + underlying_tree(t, b - 1) - underlying_tree(t, b);
                                auto w1 = (premium_tree(t, b - 1).first - calc_payoff(underlying_tree(t, b - 1))) / den;
                                auto w2 = (-premium_tree(t, b).first + calc_payoff(underlying_tree(t, b))) / den;
                                boundary[t] = w1 * underlying_tree(t, b) + w2 * underlying_tree(t, b - 1);
                            } else if (b==0) {
                                boundary[t] = underlying_tree(t,b);
                            } else {
                                //dont know, we could repeat from the next step or use nan
                                boundary[t] = boundary[t+1]*discount_factor;
                            }

                        } else if (type == instrument_type::call) {
                            for (auto it = start; it != end; it++) {
                                if (it->second) {
                                    b = it - start;
                                }
                            }
                            if(b>0 and b<t) {
                                //Not sure this is correct, the paper didnt have a formula for it.
                                auto den = premium_tree(t, b + 1).first - premium_tree(t, b).first + underlying_tree(t, b + 1) - underlying_tree(t, b);
                                auto w1 = (premium_tree(t, b + 1).first - calc_payoff(underlying_tree(t, b + 1))) / den;
                                auto w2 = (-premium_tree(t, b).first + calc_payoff(underlying_tree(t, b))) / den;
                                boundary[t] = w1 * underlying_tree(t, b) + w2 * underlying_tree(t, b + 1);
                            } else if (b==0) {
                                boundary[t] = underlying_tree(t,b);
                            } else {
                                //dont know, we could repeat from the next step or use nan
                                boundary[t] = boundary[t+1]*discount_factor;
                            }
                        }
                    }

                }

                return boundary;
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
