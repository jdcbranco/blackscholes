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
                limit /= 2;

                auto map_index = [limit] (int i) {
                    return limit-i;
                };

                auto dt = p.tau / steps;
                auto dsigma = p.sigma * sqrt(dt);
                auto u = exp(dsigma);
                auto d = exp(-dsigma);
                auto prob = (exp((p.r - p.q) * dt) - d) / (u - d);
                auto df = exp(-p.r * dt); //discount factor

                std::vector<T> und(height);
                std::vector<node<T>> layer1(height);
                std::vector<node<T>> layer2(height);

                std::vector<tf::Task> tasks;

                //This is the last layer1 (before the maturity) hence
                auto initialization_task = taskflow.for_each_index(0, height, 1, [&](int i) {
                    auto k = map_index(i);
                    T S = p.S*pow(u,k); //exp(dsigma*k);
                    und[i] = S;
                    auto p2 = p.clone(S,dt);
                    T payoff = this->instrument->payoff(S);
                    T continuation = this->blackScholes(*p2); //0.0;//
                    if(payoff > continuation) {
                        layer1[i] = std::make_pair(payoff, true);
                    } else {
                        layer1[i] = std::make_pair(continuation, false);
                    }
                });

                tasks.push_back(initialization_task);

                auto *out = &layer2;
                auto *in = &layer1;

                for (int t = steps-2; t >=0; --t) {

                    tf::Task processing_task = taskflow.for_each_index(0, height, 1, [&](int i) {
                        auto k = map_index(i);
                        T S = und[i];//p.S*pow(u,k); //this saves about 2%
                        T payoff = this->instrument->payoff(S);
                        T continuation;

                        if(i==0 or i==(height-1) or t==steps) {
                            auto p2 = p.clone(S,dt);
                            continuation = this->blackScholes(*p2);//(*in)[i].first*df; //
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

                }

                executor.run(taskflow).wait();

                std::cout << std::boolalpha;
                std::cout << "Result:" << std::endl;
                for(int i = 0; i<height; i++) {
                    std::cout << "index = " << map_index(i) << ", S = " << und[i] << ", P = " << (*in)[i].first << ", E = " << (*in)[i].second << std::endl;
                }

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
