
#include "solver.h"
#include "solver_analytical_internals.h"

#include <autodiff/forward/real.hpp>

using namespace autodiff;
using namespace bsm::internals;

namespace bsm {

    template<typename T>
    struct qdplus_method: pricing<T>, method {
        T M, N;
        qdplus_method(american const& instrument, mkt_params<double> mp): pricing<T>{instrument,mp},
        M{2*this->r/(this->sigma*this->sigma)},
        N{2*(this->r-this->q)/(this->sigma*this->sigma)}
        {}

        auto calc_price() {
            auto Sb = calculate_exercise_boundary(this->tau);
            auto S = this->S;
            auto K = this->K;
            auto tau = this->tau;
            auto r = this->r;

            if(S <= Sb) {
                return K - S;
            } else {
                auto eput = calculate_european_put<T>(this);
                auto eput_b = calculate_european_put<T>(this->clone(Sb,tau));
                auto h = 1-exp(-r*tau);
                auto qd = calc_qqd(M, N, h); //q_QD
                auto qdd = calc_qqd_deriv(M, N, h); //q_QD'(h)
                auto b = calc_b(M, N, h, qd, qdd);
                auto logSSb = log(S/Sb);
                auto c0 = calc_c0(M, N, h, qd, qdd, Sb, tau, eput);
                auto c = c0;
                return eput +
                        ((K-Sb-eput_b)/(1-b*(logSSb*logSSb)-c*logSSb))*pow(S/Sb,qd);
            }
        }

        auto get_exercise_boundary_function(T tau) {
            auto r = pricing<T>::r;
            //auto tau = pricing<T>::tau;
            auto h = 1-exp(-r*tau);
            auto qd = calc_qqd(M, N, h); //q_QD
            auto qdd = calc_qqd_deriv(M, N, h); //q_QD'(h)
            return [this,&tau,h,qd,qdd](T Sb) {
                auto& q = this->q;
                auto& K = this->K;
                auto p = this->clone(Sb,tau);
                auto d1 = calculate_d1<T>(p);
                auto eput = calculate_european_put<T>(p);
                auto c0 = calc_c0(M, N, h, qd, qdd, Sb, tau, eput);
                auto c = c0;
                return (1-exp(-q*tau)*cdf<T>(-d1))*Sb + (qd + c)*(K-Sb-eput);
            };
        }

        auto calculate_exercise_boundary(T tau) {
            T Sb = this->K;
            auto equation = get_exercise_boundary_function(tau);
            auto [u0, ux] = derivatives(equation, wrt(Sb), at(Sb));
            for(int i = 0; i<100; i++) {
                if (u0 < 1e-9) {
                    break;
                }
                Sb -= u0/ux;
            }
            return Sb;
        }

        //qdd means q_QD'(h)
        T calc_b(T M, T N, T h, T qd, T qdd) {
            return 0.5*(1-h)*M*qdd/(2*qd+N-1);
        }

        T calc_qqd(T M, T N, T h) {
            return -0.5*(N-1+sqrt((N-1)*(N-1) + 4*M/h));
        }

        T calc_qqd_deriv(T M, T N, T h) {
            return M/(h*h*sqrt((N-1)*(N-1)+ 4*M/h));
        }

        T calc_c0(T M, T N, T h, T qd, T qdd, T Sb, T tau, T eput) {
            auto& r = this->r;
            auto& K = this->K;
            return -((1-h)*M/(2*qd+N-1))*(1/h - put_theta(Sb,tau)*exp(r*tau)/(r*(K-Sb-eput)) + qdd/(2*qd+N-1) );
        }

        T put_theta(T S, T tau) {
            pricing<double> p = this->clone(S,tau);
            auto& K = this->K;
            auto& sigma = this->sigma;
            auto& r = this->r;
            auto& q = this->q;
            auto d1 = calculate_d1<double>(p);
            auto d2 = calculate_d2<double>(p);
            return r * K * exp(-r*tau) * cdf(-d2)
            - q * S * exp(-q * tau) * cdf(-d1)
            - 0.5 * sigma * S * exp(-q * tau) * pdf(d1) / sqrt(tau);
        }


    };



}