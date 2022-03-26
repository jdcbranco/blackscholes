#ifndef BSM_SOLVER_LATTICE_INTERNALS_H
#define BSM_SOLVER_LATTICE_INTERNALS_H

#include "common.h"
#include "instruments.h"
#include "solver.h"
#include "solver_analytical_internals.h"

#include <vector>
#include <iostream>
#include <algorithm>
#include <execution>
#include <utility>
#include <ranges>

#include <taskflow/taskflow.hpp>

namespace bsm {
    namespace internals {

        template<typename T>
        using node = std::pair<T,bool>;

        template<typename T>
        using calc_payoff_type = T(T const&);

        template<typename T>
        struct superpositioned_binomial_lattice_method: american_method {
            const int steps;
            std::shared_ptr<american> instrument;
            std::function<pricing_function<T>> blackScholes;

            superpositioned_binomial_lattice_method(american_put const& instrument,int steps):
                steps{steps}, instrument{std::make_shared<american_put>(instrument)} {
                blackScholes = [](pricing<T> const& p) { return calculate_european_put<T>(p); };
            }
            superpositioned_binomial_lattice_method(american_call const& instrument,int steps):
                    steps{steps}, instrument{std::make_shared<american_call>(instrument)} {
                blackScholes = [](pricing<T> const& p) { return calculate_european_call<T>(p); };
            }
            superpositioned_binomial_lattice_method(superpositioned_binomial_lattice_method const&) = default;
            superpositioned_binomial_lattice_method(superpositioned_binomial_lattice_method &&) noexcept = default;


            void solve(pricing<T> const& p) {

                tf::Executor executor;
                tf::Taskflow taskflow;

                //This controls the height of the lattice
//                auto max_S = p.S * exp(p.sigma * sqrt(p.tau) * 4);
//                auto min_S = p.S * exp(p.sigma * sqrt(p.tau) * -8);
//                std::cout << "Max S = " << max_S << ", Min S = " << min_S << std::endl;

                //int limit_up = 6*sqrt(steps);
                int limit = 8*sqrt(steps);
                int height = 2 * limit + 1;
                if(instrument->type==instrument_type::put)
                    limit /= 2;
                else if(instrument->type==instrument_type::call)
                    limit *= 1.5;

                auto map_index = [limit] (int i) {
                    return limit-i;
                };

                auto dt = p.tau / steps;
                auto dsigma = p.sigma * sqrt(dt);
                auto u = exp(dsigma);
                auto d = exp(-dsigma);
                auto prob = (exp((p.r - p.q) * dt) - d) / (u - d);
                auto df = exp(-p.r * dt); //discount factor

                std::vector<T> underlying(height);
                std::vector<node<T>> layer1(height);
                std::vector<node<T>> layer2(height);
                std::vector<T> boundary(steps);

                std::vector<tf::Task> tasks;

                //This is the last layer1 (before the maturity) hence
                auto initialization_task = taskflow.for_each_index(0, height, 1, [&](int i) {
                    auto k = map_index(i);
                    T S = p.S*pow(u,k); //exp(dsigma*k);
                    underlying[i] = S;
                    auto p2 = p.clone(S,dt);
                    T payoff = this->instrument->payoff(S);
                    T continuation = this->blackScholes(*p2);
                    if(payoff > continuation) {
                        layer1[i] = std::make_pair(payoff, true);
                    } else {
                        layer1[i] = std::make_pair(continuation, false);
                    }
                });

                tasks.push_back(initialization_task);

                auto *out = &layer2;
                auto *in = &layer1;

                int tval = steps-2;
                int *tptr = &tval;

                for (int t = steps-2; t>=0; --t) {

                    tf::Task processing_task = taskflow.for_each_index(0, height, 1, [&](int i) {
                        auto k = map_index(i);
                        T S = underlying[i];//p.S*pow(u,k); //this saves about 2%
                        T payoff = this->instrument->payoff(S);
                        T continuation;

                        if(i==0 or i==(height-1)) {
                            auto p2 = p.clone(S,dt);
                            continuation = this->blackScholes(*p2);
                        } else {
                            T premium_up = (*in)[i-1].first;
                            T premium_down = (*in)[i+1].first;
                            continuation = (prob*premium_up + (1.0-prob)*premium_down) * df;
                        }

                        if(payoff > continuation) {
                            (*out)[i] = std::make_pair(payoff, true);
                        } else {
                            (*out)[i] = std::make_pair(continuation,false);
                        }

                    });

                    tasks.back().precede(processing_task);
                    tasks.push_back(processing_task);

                    tf::Task swap_task = taskflow.emplace([&](){
                        std::swap(in,out);
                    });

                    tasks.back().precede(swap_task);
                    tasks.push_back(swap_task);

                    tf::Task boundary_task = taskflow.emplace([&]() {
                        int b = -1;
                        int t = *tptr;
                        if(instrument->type==instrument_type::put) {
                            auto B = std::max_element(std::execution::par_unseq, out->begin(), out->end(),
                                [](node<T> const &a, node<T> const &b) {
                                  if (a.second == false and b.second == true) {
                                      return true;
                                  } else if (a.second == false and b.second == false) {
                                      return a.first > b.first;
                                  } else if (a.second == true and b.second == false) {
                                      return false;
                                  } else if (a.second == true and b.second == true) {
                                      return a.first > b.first;
                                  }
                                });
                            b = B - out->begin();
                            //b = 10;
                            if(b>0) {
                                auto& premium = *in;
                                //This approximation is based on paper "Discrete and continuous time approximations of the optiomal exercise boundary of American options - Basso, Nardon, Pianca"
                                auto den = premium[b-1].first - premium[b].first + underlying[b-1] - underlying[b];
                                auto w1 = (premium[b-1].first - this->instrument->payoff(underlying[b-1])) / den;
                                auto w2 = (-premium[b].first + this->instrument->payoff(underlying[b])) / den;
                                boundary[t] = w1 * underlying[b] + w2 * underlying[b-1];
                            } else if (b==0) {
                                boundary[t] = underlying[b];
                            } else {
                                //don't know, we could repeat from the next step
                                boundary[t] = boundary[t+1]*df;
                            }
                        } else if(instrument->type==instrument_type::call) {
                            auto B = std::max_element(std::execution::par_unseq, out->begin(), out->end(),
                                [](node<T> const &a, node<T> const &b) {
                                  if (a.second == false and b.second == true) {
                                      return true;
                                  } else if (a.second == false and b.second == false) {
                                      return a.first < b.first;
                                  } else if (a.second == true and b.second == false) {
                                      return false;
                                  } else if (a.second == true and b.second == true) {
                                      return a.first < b.first;
                                  }
                                });
                            b = B - out->begin();
                            if(b>0 and b<t) {
                                //Not sure this is correct, the paper didnt have a formula for it.
                                auto den = (*in)[b+1].first - (*in)[b].first + underlying[b+1] - underlying[b];
                                auto w1 = ((*in)[b+1].first - this->instrument->payoff(underlying[b+1])) / den;
                                auto w2 = (-(*in)[b].first + this->instrument->payoff(underlying[b])) / den;
                                boundary[t] = w1 * underlying[b] + w2 * underlying[b+1];
                            } else if (b==0) {
                                boundary[t] = underlying[b];
                            } else {
                                //don't know, we could repeat from the next step or use nan
                                boundary[t] = boundary[t+1]*df;
                            }
                        }
                        (*tptr) -= 1;
                    });

                    tasks.back().precede(boundary_task);
                    tasks.push_back(boundary_task);

                }

                executor.run(taskflow).wait();

                std::cout << std::boolalpha;
                std::cout << "Result:" << std::endl;
                for(int i=0; i<height; i++) {
                    std::cout << "index = " << map_index(i) << ", S = " << underlying[i] << ", P = " << (*in)[i].first << ", E = " << (*in)[i].second << std::endl;
                }
                std::cout << "Boundary:" << std::endl;
                for(int i=0; i<steps; i++) {
                    std::cout << "Sb["<< i <<"] = " << boundary[i] << std::endl;
                }

                //TODO Complete this

            }

            double price() override {
                return 0;
            }

            double delta() override {
                return 0;
            }

            double gamma() override {
                return 0;
            }

            double vega() override {
                return 0;
            }

            double theta() override {
                return 0;
            }

            double rho() override {
                return 0;
            }

            double psi() override {
                return 0;
            }

            long double exercise_boundary(long double _tau) override {
                return 0;
            }

        };

    }
}

#endif //BSM_SOLVER_LATTICE_INTERNALS_H
