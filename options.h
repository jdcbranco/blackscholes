#ifndef BLACKSCHOLES_OPTIONS_H
#define BLACKSCHOLES_OPTIONS_H
#include <algorithm>
#include <ranges>
#include <execution>
#include <vector>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <utility>

namespace options {
    using namespace std;

    template<typename T>
    struct params {
        T S;
        T sigma;
        T tau;
        T r;
        T q;
    public:
        params(T S, T sigma, T tau, T r, T q): S{S}, sigma{sigma}, tau{tau}, r{r}, q{q}
        {
        }

        template<std::size_t Index>
        std::tuple_element_t<Index, params<T>>& get() &
        {
            if constexpr (Index == 0) return S;
            if constexpr (Index == 1) return sigma;
            if constexpr (Index == 2) return tau;
            if constexpr (Index == 3) return r;
            if constexpr (Index == 4) return q;
        }

        template<std::size_t Index>
        std::tuple_element_t<Index, params<T>> const& get() const&
        {
            if constexpr (Index == 0) return S;
            if constexpr (Index == 1) return sigma;
            if constexpr (Index == 2) return tau;
            if constexpr (Index == 3) return r;
            if constexpr (Index == 4) return q;
        }

        template<std::size_t Index>
        std::tuple_element_t<Index, params<T>>& get() &&
        {
            if constexpr (Index == 0) return move(S);
            if constexpr (Index == 1) return move(sigma);
            if constexpr (Index == 2) return move(tau);
            if constexpr (Index == 3) return move(r);
            if constexpr (Index == 4) return move(q);
        }

        template<std::size_t Index>
        std::tuple_element_t<Index, params<T>> const& get() const&&
        {
            if constexpr (Index == 0) return move(S);
            if constexpr (Index == 1) return move(sigma);
            if constexpr (Index == 2) return move(tau);
            if constexpr (Index == 3) return move(r);
            if constexpr (Index == 4) return move(q);
        }

    };

    template<typename T>
    using market_variables = shared_ptr<options::params<T>>;

    template<typename I, typename T>
    concept Exercisable = requires(I instrument, T price_at_maturity) {
        { instrument.payoff(price_at_maturity) } -> std::convertible_to<T>;
    };

    template<typename I, typename T>
    concept ClosedFormSolution = requires(I instrument, params<T> const& params) {
        { instrument.closed_form_price(params) } -> std::convertible_to<T>;
    };

    /**
     * Standard Normal CDF
     * @param x
     * @return
     */
    template<typename T>
    T phi(T const& x);

    template<typename T = double>
    class forward {
        T K;
    public:
        forward(T strike): K{strike} {}
        T payoff(T price_at_maturity) const {
            return price_at_maturity - K;
        }
        T closed_form_price(params<T> const& params) {
            auto& [S, sigma, tau, r, q] = params;
            return S * exp(-q*tau) -  K * exp(- r*tau);
        }
        T fair_price(params<T> const& params) {
            auto& [S, sigma, tau, r, q] = params;
            return S * exp ((r-q) * tau);
        }
    };

    template<typename T = double>
    class european {
    public:
        european(T strike) : strike_{strike} {}
        virtual T payoff(T price_at_maturity) const = 0;
        T strike() const { return strike_; }
    private:
        T strike_;
    };

    template<typename T>
    class european_call: public european<T> {
    public:
        european_call(T strike): european<T>{strike} {}
        T payoff(T price_at_maturity) const override {
            return max(price_at_maturity - this->strike(),0.0);
        }
        T closed_form_price(params<T> const& params) {
            auto& [S, sigma, tau, r, q] = params;
            auto K = european<T>::strike();
            T const d1 = (log(S / K) + (r - q + sigma*sigma / 2) * tau) / (sigma * sqrt(tau));
            T const d2 = (log(S / K) + (r - q - sigma*sigma / 2) * tau) / (sigma * sqrt(tau));
            return exp(-q * tau) * S * phi(d1) - exp(-r * tau) * K * phi(d2);
        }
    };

    template<typename T>
    class european_put: public european<T> {
    public:
        european_put(T strike): european<T>{strike} {}
        T payoff(T price_at_maturity) const override {
            return max(this->strike()-price_at_maturity,0.0);
        }
        T closed_form_price(params<T> const& params) {
            auto& [S, sigma, tau, r, q] = params;
            auto K = european<T>::strike();
            T const d1 = (log(S / K) + (r - q + sigma*sigma / 2) * tau) / (sigma * sqrt(tau));
            T const d2 = (log(S / K) + (r - q - sigma*sigma / 2) * tau) / (sigma * sqrt(tau));
            return exp(-r * tau) * K * phi(-d2) - exp(-q * tau) * S * phi(-d1);
        }
    };

    template<typename T = double>
    class american {
        T K;
    public:
        american(T strike): K{strike} {}
        virtual T payoff(T price) const = 0;
        T strike() const { return K; }
    };

    template<typename T>
    class american_call: public american<T> {
    public:
        american_call(T strike): american<T>{strike} {}
        T payoff(T price) const override {
            return max(price - this->strike(), 0.0);
        }
    };

    template<typename T>
    class american_put: public american<T> {
    public:
        american_put(T strike): american<T>{strike} {}
        T payoff(T price) const override {
            return max(this->strike()-price, 0.0);
        }
    };

    template<typename T>
    class pricing {
        T price_;
        market_variables<T> params_ptr;
    public:
        pricing(T price, market_variables<T> const& params):
                price_{price},
                params_ptr{params} { }
        pricing(T && price, market_variables<T> const& params):
                price_{price},
                params_ptr{params} { }
        pricing(pricing<T> const&) = default;
        pricing(pricing<T> &&) noexcept = default;

        explicit operator T const& () { return price_; }

        inline auto price() const { return price_; }
        inline auto& variables() const { return *params_ptr; }
        inline auto& S() const { return params_ptr->S; }
        inline auto& sigma() const { return params_ptr->sigma; }
        inline auto& tau() const { return params_ptr->tau; }
        inline auto& r() const { return params_ptr->r; }
        inline auto& q() const { return params_ptr->q; }
    };


    template<typename T>
    ostream& operator<<(ostream& out, pricing<T> pricing) {
        out << (T)pricing;
        return out;
    }



}

namespace std {
    template<typename T>
    struct tuple_size<options::params<T>> : std::integral_constant<size_t, 5> { };

    template<size_t Index, typename T>
    struct tuple_element<Index, options::params<T>>
            : tuple_element<Index, tuple<T, T, T, T, T>>
    {
    };
}


#endif //BLACKSCHOLES_OPTIONS_H
