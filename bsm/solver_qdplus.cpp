#include "solver.h"
#include "solver_analytical_internals.h"

#include <autodiff/forward/dual.hpp>

using namespace autodiff;
using namespace bsm::internals;

namespace bsm {

    template<typename T>
    using exercise_boundary_function = T(T const&);

    //This is still a working in progress. The code below needs a lot of work, but has been debugged and I know it works.
    //Paper: Analytical Approximations for the Critical Stock Prices of American Options: A Performance Comparison
    //by Minqiang Li, Li (2009)
    //This model (QD+) is considered to provide good approximation for the exercise boundary.

    template<typename T>
    struct qdplus_method_core: pricing<T> {
        T M, N;
        qdplus_method_core(american const& instrument, mkt_params<double> mp): pricing<T>{instrument,mp},
        M{2.0*this->r/(this->sigma*this->sigma)},
        N{2.0*(this->r-this->q)/(this->sigma*this->sigma)}
        {}

        T calc_price() {
            auto Sb = calculate_exercise_boundary(this->tau);
            auto S = this->S;
            auto K = this->K;
            auto tau = this->tau;
            auto r = this->r;

            if(S <= Sb) {
                return K - S;
            } else {
                auto eput = calculate_european_put<T>(*this);
                auto eput_b = calculate_european_put<T>(*(this->clone(Sb,tau)));
                T h = 1.0-exp(-r*tau);
                auto qd = calc_qqd(M, N, h); //q_QD
                auto qdd = calc_qqd_deriv(M, N, h); //q_QD'(h)
                auto b = calc_b(M, N, h, qd, qdd);
                auto logSSb = log(S/Sb);
                auto c0 = calc_c0(M, N, h, qd, qdd, Sb, tau, eput_b);
                auto c = c0;
                return eput +
                        ((K-Sb-eput_b)/(1.0-b*(logSSb*logSSb)-c*logSSb))*pow(S/Sb,qd);
            }
        }

        std::function<exercise_boundary_function<T>> get_exercise_boundary_function(T tau) {
            auto r = pricing<T>::r;
            //auto tau = pricing<T>::tau;
            T h = 1.0-exp(-r*tau);
            auto qd = calc_qqd(M, N, h); //q_QD
            auto qdd = calc_qqd_deriv(M, N, h); //q_QD'(h)
            return [this,tau,h,qd,qdd](T const& Sb) {
                auto q = this->q;
                auto K = this->K;
                auto p = this->clone(Sb,tau);
                auto d1 = calculate_d1<T>(*p);
                auto eput_b = calculate_european_put<T>(*p);
                auto c0 = calc_c0(M, N, h, qd, qdd, Sb, tau, eput_b);
                auto c = c0;
                return (1.0-exp(-q*tau)*cdf<T>(-d1))*Sb + (qd + c)*(K - Sb - eput_b);
            };
        }

        T calculate_exercise_boundary(T tau) {
            T Sb = this->K;
            auto equation = get_exercise_boundary_function(tau);

            for(int i = 0; i<100; i++) {
                auto [u0, ux] = derivatives(equation, wrt(Sb), at(Sb));
                if (u0 < 1e-9) {
                    break;
                }
                Sb -= u0/ux;//val(u0/ux);
            }
            return Sb;
        }

        //qdd means q_QD'(h)
        T calc_b(T M, T N, T h, T qd, T qdd) {
            return 0.5*(1.0-h)*M*qdd/(2.0*qd+N-1.0);
        }

        T calc_qqd(T M, T N, T h) {
            return -0.5*(N-1+sqrt((N-1.0)*(N-1.0) + 4.0*M/h));
        }

        T calc_qqd_deriv(T M, T N, T h) {
            return M/(h*h*sqrt((N-1.0)*(N-1.0)+ 4*M/h));
        }

        T calc_c0(T M, T N, T h, T qd, T qdd, T Sb, T tau, T eput) {
            auto r = this->r;
            auto K = this->K;
            return -((1.0-h)*M/(2.0*qd+N-1.0))*(1.0/h - put_theta(Sb,tau)*exp(r*tau)/(r*(K-Sb-eput)) + qdd/(2.0*qd+N-1.0) );
        }

        T put_theta(T S, T tau) {
            auto p = this->clone(S,tau);
            auto& K = this->K;
            auto& sigma = this->sigma;
            auto& r = this->r;
            auto& q = this->q;
            auto d1 = calculate_d1<T>(*p);
            auto d2 = calculate_d2<T>(*p);
            return r * K * exp(-r*tau) * cdf<T>(-d2)
            - q * S * exp(-q * tau) * cdf<T>(-d1)
            - 0.5 * sigma * S * exp(-q * tau) * pdf<T>(d1) / sqrt(tau);
        }

    };

    struct qdplus_method: pricing<double>, american_method {
        qdplus_method_core<dual> core;
        qdplus_method(american const& instrument, mkt_params<double> mp): pricing{instrument,mp}, core{instrument,mp} {

        }

        double price() override {
            return val(core.calc_price());
        }

        double delta() override {
            return NAN;
        }

        double gamma() override {
            return NAN;
        }

        double vega() override {
            return NAN;
        }

        double theta() override {
            return NAN;
        }

        double rho() override {
            return NAN;
        }

        double psi() override {
            return NAN;
        }

        double exercise_boundary(double _tau) override {
            return val(core.calculate_exercise_boundary(_tau));
        }
    };

    template<>
    std::unique_ptr<american_method> qdplus_solver<autodiff_off>::operator()(american_put& instrument) {
        qdplus_method gp{instrument, mktParams};
        return std::make_unique<qdplus_method>(gp);
    }

}