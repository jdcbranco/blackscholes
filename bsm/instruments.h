#ifndef BSM_INSTRUMENTS_H
#define BSM_INSTRUMENTS_H

#include "common.h"
#include "chrono.h"

#include <concepts>
#include <chrono>
#include <iostream>

namespace bsm {
    using namespace bsm::chrono;

    struct instrument {
        const long double K;
        const datetime maturity;
        instrument(long double const& strike, datetime const& maturity): K{strike}, maturity{maturity} {}
    };

    template<typename I>
    concept Expirable = requires(I instrument) {
        {  instrument.maturity } -> std::convertible_to<datetime>;
    };

    template<typename I, typename T>
    concept Exercisable = requires(I instrument, T price_at_maturity) {
        { instrument.payoff(price_at_maturity) } -> std::convertible_to<T>;
    };

    template<typename I, typename T>
    concept ClosedFormSolution = requires(I instrument, mkt_params<T> const& params) {
        { instrument.closed_form_price(params) } -> std::convertible_to<T>;
    };

    struct european: instrument {
        inline european(long double strike, datetime const& maturity) : instrument{strike, maturity} {}
        virtual long double payoff(long double price_at_maturity) const = 0;
    };

    struct european_forward: european {
        inline european_forward(long double const& strike, datetime const& maturity): european{strike, maturity} {}
        inline european_forward(european_forward const& other) = default;
        inline european_forward(european_forward && other) noexcept = default;
        inline long double payoff(long double price_at_maurity) const override {
            return price_at_maurity - K;
        }
    };

    struct european_call: european {
        inline european_call(long double strike, datetime const& maturity): european{strike, maturity} {}
        inline european_call(long double strike, year_month_day const& maturity): european{strike, maturity} {}
        inline long double payoff(long double price_at_maturity) const override {
            return std::max(price_at_maturity-K,0.0L);
        }
    };

    struct european_put: european {
        inline european_put(long double strike, datetime const& maturity): european{strike, maturity} {}
        inline european_put(long double strike, year_month_day const& maturity): european{strike, maturity} {}
        inline long double payoff(long double price_at_maturity) const override {
            return std::max(K-price_at_maturity,0.0L);
        }
    };

    struct american: instrument {
        inline american(long double strike, datetime const& maturity): instrument{strike,maturity} {} // K{strike}, maturity{maturity} {}
        virtual long double payoff(long double price) const = 0;
    };

    struct american_call: american {
        american_call(long double strike, datetime const& maturity): american{strike, maturity} {}
        american_call(long double strike, year_month_day const& maturity): american{strike, maturity} {}
        inline long double payoff(long double price) const override {
            return std::max(price - K, 0.0L);
        }
    };

    struct american_put: american {
        american_put(long double strike, datetime const& maturity): american{strike, maturity} {}
        american_put(long double strike, year_month_day const& maturity): american{strike, maturity} {}
        inline long double payoff(long double price) const override {
            return std::max(K - price, 0.0L);
        }
    };

}

#endif //BSM_INSTRUMENTS_H
