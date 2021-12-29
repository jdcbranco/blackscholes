#ifndef BLACKSCHOLES_CRR_H
#define BLACKSCHOLES_CRR_H

#include <algorithm>
#include <ranges>
#include <execution>
#include <vector>
#include <cassert>
#include <iostream>
#include <iomanip>

#include "lattice.h"
#include "options.h"

namespace crr {
    using namespace std;
    using namespace options;
    using namespace lattice;

    template<typename T>
    class crrmethod {
        bintree<T> underlying_tree;
        bintree<T> premium_tree;
        //options::params<T> params_;
        market_variables<T> params_ptr;
        T u_, d_, p_, discount_factor_;
    public:
        crrmethod(params<T> parameters, int steps): params_ptr{std::make_shared<params<T>>(parameters)},
        underlying_tree{steps + 1}, premium_tree{steps + 1} {
            auto dt = params_ptr->tau / steps;
            auto sqrt_dt = sqrt(dt);
            //auto sigma = sqrt(params.variance);
            auto u = u_ = exp(params_ptr->sigma * sqrt_dt);
            auto d = d_ = exp(-params_ptr->sigma * sqrt_dt);
            p_ = (exp((params_ptr->r - params_ptr->q)*dt) - d) / (u - d);
            cout << "CRR tree up prob = " << p_ << endl;
            discount_factor_ = exp(-params_ptr->r * dt);
            auto S = params_ptr->S;
            underlying_tree.set(0, 0, S);//TODO Implement some syntatic sugar: underlying_tree(0,0) = S;
            vector<int> indices(steps+1);
            iota(indices.begin(),indices.end(), 0);
            for (int t = 1; t <= steps; ++t) {
                auto start = indices.begin();
                auto end = start+t+1;
                auto [output, _ ] = underlying_tree(t);
                transform(execution::par_unseq, start, end, output, [&S,&u,&d,&t](int i) {
                    return S*pow(u,t-i)*pow(d,i);
                });
            }
        }

        template<typename I> requires options::Exercisable<I,T>
        auto solve(I instrument) {
            auto steps = premium_tree.size();
            auto last_t = steps-1;
            auto p = p_;
            auto discount_factor = discount_factor_;

            {
                auto [start, end] = underlying_tree(last_t);
                auto [output, _] = premium_tree(last_t);
                transform(execution::par_unseq, start, end, output,
                          [&instrument](T S_t) {
                              return instrument.payoff(S_t);
                          });
            }

            vector<int> indices(steps);
            iota(indices.begin(),indices.end(), 0);
            for(int t = last_t-1; t>=0; t--) {
                auto start = indices.begin();
                auto end = start+t+1;

                auto [output, _] = premium_tree(t);
                auto premium_next_step = get<0>(premium_tree(t+1));

                transform(execution::par_unseq, start, end, output,
                          [premium_next_step, p, discount_factor](int i) {
                    auto premium_up = *(premium_next_step+i);
                    auto premium_down = *(premium_next_step+i+1);
                    return (p*premium_up + (1.0-p)*premium_down)*discount_factor;
                });
            }

            auto price = premium_tree.root();
            options::pricing<T> pricing{price, params_ptr};
            return pricing;
        }

        auto underlying() const {
            return underlying_tree;
        }

        auto premium() const {
            return premium_tree;
        }

    };

    template<typename T>
    ostream& operator<<(ostream& out, crrmethod<T> const& crrtree) {
        out << crrtree.underlying();
        return out;
    }

}

#endif //BLACKSCHOLES_CRR_H
