#ifndef BLACKSCHOLES_LSMC_H
#define BLACKSCHOLES_LSMC_H

#include <execution>
#include <Eigen/Dense>
#include "options.h"
#include "random.h"

namespace lsmc {
    using namespace options;

    template<typename T = double>
    class lsmc_method {
    private:
        market_variables <T> params_ptr;
        random_normal <double> random;
    public:
        explicit inline lsmc_method(params <T> const &parameters) : params_ptr{std::make_shared<params<T>> (parameters)} {}

        template<typename I> requires Exercisable<I,T>
        pricing<T> solve(I instrument, int steps, int paths) {
            auto &[S, sigma, tau, r, q] = (*params_ptr);

            auto dt = tau / steps;
            auto sqrt_dt = sqrt(dt);
            auto R = exp((r-q)*dt);
            auto variance = sigma*sigma;

            Eigen::MatrixX<T> X(paths,steps+1);

            for(int p=0; p<paths; p++) {
                T x = S;
                auto z = random(steps+1);
                for(int s=0; s<steps+1; s++) {
                    x *= R * exp(-variance*dt/2.0 + sigma*sqrt_dt*z[s]);
                    X(p,s) = (double)x;
                }
            }

            Eigen::VectorX<T> V(paths);
            for(int p=0; p<paths; p++) {
                V[p] = instrument.payoff(X(p,steps));
            }

            //backwards
            std::vector<int> itm;
            std::vector<bool> is_itm;
            for(int s=steps-1; s>= 0; s--) {
                //std::cout << "Step " << s << "\n";

                itm.clear();
                is_itm.clear();
                //identify where the Vs are in the money
                for(int p=0; p<paths; p++) {
                    if(V[p] > 0.0) {
                        itm.push_back(p);
                        is_itm.push_back(true);
                    } else{
                        is_itm.push_back(false);
                    }
                }

                Eigen::VectorX<T> b(itm.size());
                Eigen::MatrixX<T> A(itm.size(),3);

                A(Eigen::all,0).setOnes();
                for(int i = 0 ; i<itm.size(); i++) {
                    b(i) = V[itm[i]];
                    A(i,1) = X(itm[i], s);
                    A(i,2) = A(i,1)*A(i,1);
                }

                Eigen::VectorX<T> coefficients = A.colPivHouseholderQr().solve(b);
                Eigen::VectorX<T> continuation = A * coefficients;

                for(int p=0,itm_idx=0; p<paths; p++) {
                    auto exercise_value = instrument.payoff(X(p,s)) / R;
                    //if itm and exercise value > continuation
                    if(is_itm[p] && exercise_value > continuation[itm_idx++]) {
                        V[p] = exercise_value / R;
                    } else {
                        V[p] = V[p] / R;
                    }
                }
            }

            auto price = max(V.mean(), instrument.payoff(S));
            pricing<T> pricing{price, params_ptr};
            return pricing;
        }
    };

}

#endif //BLACKSCHOLES_LSMC_H
