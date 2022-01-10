#include "solver.h"
#include "solver_analytical_internals.h"
#include "solver_american_internals.h"

#include <autodiff/forward/dual.hpp>

using namespace autodiff;
using namespace bsm::internals;

namespace bsm {

    template<typename T>
    using exercise_boundary_function = T(T const&);

    //Experimenting wiht long double duals
    using ldual = HigherOrderDual<3, long double>;

    /**
     * Paper: Analytical Approximations for the Critical Stock Prices of American Options: A Performance Comparison
       by Minqiang Li, Li (2009). This QD+ model is considered to provide good approximation for the exercise boundary.
     * @tparam T
     */
    template<typename T>
    struct qdplus_method_core: pricing<T> {
        T M, N;
        bool call;

        qdplus_method_core(pricing<T> const& p, bool call = false):
                pricing<T>{p.S,p.K,p.sigma,p.tau,p.r,p.q},
                call{call},
                M{2.0L*this->r/(this->sigma*this->sigma)},
                N{2.0L*(this->r-this->q)/(this->sigma*this->sigma)}
        {}

        /**
         * Calculates price based on exercise boundary
         * @param Sb Exercise boundary
         * @return
         */
        T calc_price(T const& Sb) {
            auto& S = this->S;
            auto& K = this->K;
            auto& tau = this->tau;
            auto& r = this->r;

            if(call and S >= Sb) {
                return S - K;
            } else if(not call and S <= Sb) {
                return K - S;
            } else {
                auto european_option = call ? calculate_european_call<T>(*this) : calculate_european_put<T>(*this);
                if (never_optimal_exercise<T>(*this,call))
                    return european_option;
                else {
                    auto european_option_at_boundary = call ? calculate_european_call<T>(*(this->clone(val(Sb), tau))) : calculate_european_put<T>(*(this->clone(val(Sb), tau)));
                    T h = 1.0L - exp(-r * tau);
                    auto qd = calc_qqd(M, N, h); //q_QD
                    auto qdd = calc_qqd_deriv(M, N, h); //q_QD'(h)
                    auto b = calc_b(M, N, h, qd, qdd);
                    auto logSSb = log(S / Sb);
                    auto c0 = calc_c0(M, N, h, qd, qdd, Sb, tau, european_option_at_boundary);
                    auto c = c0;
                    return european_option +
                           ((K - Sb - european_option_at_boundary) / (1.0 - b * (logSSb * logSSb) - c * logSSb)) * pow(S / Sb, qd);
                }
            }
        }

        std::function<exercise_boundary_function<T>> get_exercise_boundary_function(T const& tau) {
            auto& r = this->r;
            auto h = 1.0L-exp(-r*tau);
            auto qd = calc_qqd(M, N, h); //q_QD
            auto qdd = calc_qqd_deriv(M, N, h); //q_QD'(h)
            if(call)
                return [this,tau,h,qd,qdd](T const& Sb) {
                    auto& q = this->q;
                    auto& K = this->K;
                    auto p = this->clone(Sb,tau);
                    auto d1 = calculate_d1<T>(*p);
                    auto ecall_b = calculate_european_call<T>(*p);
                    auto c0 = calc_c0(M, N, h, qd, qdd, Sb, tau, ecall_b);
                    auto c = c0;
                    return abs((1.0-exp(-q*tau)*cdf<T>(d1))*Sb - (qd + c)*(Sb - K - ecall_b));
                };
            else
                return [this,tau,h,qd,qdd](T const& Sb) {
                    auto& q = this->q;
                    auto& K = this->K;
                    auto p = this->clone(Sb,tau);
                    auto d1 = calculate_d1<T>(*p);
                    auto eput_b = calculate_european_put<T>(*p);
                    auto c0 = calc_c0(M, N, h, qd, qdd, Sb, tau, eput_b);
                    auto c = c0;
                    return abs((1.0-exp(-q*tau)*cdf<T>(-d1))*Sb + (qd + c)*(K - Sb - eput_b));
                };
        }

        T calculate_exercise_boundary(T const& tau) {
            auto& r = this->r;
            auto& q = this->q;
            if(never_optimal_exercise<T>(*this,call)) {
                return call? INFINITY : 0.0;
            }
            if(tau==0) {
                return exercise_boundary_at_maturity<T>(*this,call ? instrument_type::call : instrument_type::put);
            }

            T Sb = this->K;
            auto equation = get_exercise_boundary_function(tau);

            for(int i = 0; i<100; i++) {
                auto [u0, ux, uxx, uxxx] = derivatives(equation, wrt(Sb), at(Sb));
                if (u0 < 1e-9) {
                    break;
                }
                Sb -= u0/ux;
            }
            std::cout << (call?"Call ":"Put ") << "Sb = " << Sb << ", K = " << this->K << std::endl;
            return Sb;
        }

        //qdd means q_QD'(h)
        T calc_b(T const& M, T const& N, T const& h, T const& qd, T const& qdd) {
            return 0.5*(1.0-h)*M*qdd/(2.0*qd+N-1.0);
        }

        T calc_qqd(T const& M, T const& N, T const& h) {
            if(call) {
                return -0.5L*(N-1.0L-sqrt((N-1.0L)*(N-1.0L) + 4.0L*M/h));
            } else {
                return -0.5L*(N-1.0L+sqrt((N-1.0L)*(N-1.0L) + 4.0L*M/h));
            }
        }

        T calc_qqd_deriv(T const& M, T const& N, T const& h) {
            return M/(h*h*sqrt((N-1.0L)*(N-1.0L)+ 4.0L*M/h));
        }

        T calc_c0(T const& M, T const& N, T const& h, T const& qd, T const& qdd, T const& Sb, T const& tau, T const& eput_b) {
            auto& r = this->r;
            auto& K = this->K;
            return -((1.0-h)*M/(2.0*qd+N-1.0))*(1.0/h - put_theta(Sb,tau)*exp(r*tau)/(r*(K - Sb - eput_b)) + qdd / (2.0 * qd + N - 1.0) );
        }

        T put_theta(T const& S, T const& tau) {
            auto p = this->clone(S,tau);
            return calculate_theta<T>(*p,-1.0);
        }

    };

    struct qdplus_method: american_method {
        protected:
        bool call;
        ldual price_;
        ldual delta_;
        ldual gamma_;
        ldual vega_;
        ldual theta_;
        ldual rho_;
        ldual psi_;
        public:
            pricing<ldual> dp;

            qdplus_method(american_call const& instrument, mkt_params<double> mp, bool symmetric = false): dp{instrument,mp}, call{true} {
                //All greeks are calculated wrt to this exercise boundary.
                long double Sb = exercise_boundary(val(dp.tau));

                auto calc = [Sb](pricing<ldual> const& p) {
                    qdplus_method_core<ldual> core_{p,true};
                    return core_.calc_price(Sb);
                };

                auto [_price, _delta, _gamma, _] = symmetric ? derivatives(calc, wrt(dp.K, dp.K), at(dp)) : derivatives(calc, wrt(dp.S, dp.S), at(dp));
                price_ = _price;
                delta_ = _delta;
                gamma_ = _gamma;
                vega_ = derivative(calc, wrt(dp.sigma), at(dp));
                theta_ = derivative(calc, wrt(dp.tau), at(dp));
                rho_ = derivative(calc, wrt(symmetric? dp.q: dp.r), at(dp));
                psi_ = derivative(calc, wrt(symmetric? dp.r: dp.q), at(dp));
            }

            qdplus_method(american_put const& instrument, mkt_params<double> mp, bool symmetric = false): dp{instrument,mp}, call{false} {
                //All greeks are calculated wrt to this exercise boundary.
                long double Sb = exercise_boundary(val(dp.tau));

                auto calc = [Sb](pricing<ldual> const& p) {
                    qdplus_method_core<ldual> core_{p,false};
                    return core_.calc_price(Sb);
                };

                auto [_price, _delta, _gamma, _] = symmetric ? derivatives(calc, wrt(dp.K, dp.K), at(dp)) : derivatives(calc, wrt(dp.S, dp.S), at(dp));
                price_ = _price;
                delta_ = _delta;
                gamma_ = _gamma;
                vega_ = derivative(calc, wrt(dp.sigma), at(dp));
                theta_ = derivative(calc, wrt(dp.tau), at(dp));
                rho_ = derivative(calc, wrt(symmetric? dp.q: dp.r), at(dp));
                psi_ = derivative(calc, wrt(symmetric? dp.r: dp.q), at(dp));
            }

            double price() override {
                return val(price_);
            }

            double delta() override {
                return val(delta_);
            }

            double gamma() override {
                return val(gamma_);
            }

            double vega() override {
                return val(vega_);
            }

            double theta() override {
                return -val(theta_);
            }

            double rho() override {
                return val(rho_);
            }

            double psi() override {
                return val(psi_);
            }

            long double exercise_boundary(long double _tau) override {
                qdplus_method_core<ldual> core_{dp,call};
                return val(core_.calculate_exercise_boundary(_tau));
            }
    };

    template<>
    std::unique_ptr<american_method> qdplus_solver<autodiff_off>::operator()(american_put& instrument) {
        qdplus_method gp{instrument, mktParams};
        return std::make_unique<qdplus_method>(gp);
    }

    template<>
    std::unique_ptr<american_method> qdplus_solver<autodiff_off>::operator()(american_call& instrument) {
        qdplus_method gp{instrument, mktParams};
        return std::make_unique<qdplus_method>(gp);
    }

}