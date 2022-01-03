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
        const double K;
        const datetime maturity;
        instrument(double const& strike, datetime const& maturity): K{strike}, maturity{maturity} {}
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

    struct forward: instrument {
        inline forward(double const& strike, datetime const& maturity): instrument{strike,maturity} {} //K{strike}, maturity{maturity} {}
        forward(forward const& other) = default;
        forward(forward && other) noexcept = default;
        inline double payoff(double price_at_maurity) const {
            return price_at_maurity - K;
        }
    };

    struct european: instrument {
        inline european(double strike, datetime const& maturity) : instrument{strike, maturity} {} //K{strike}, maturity{maturity} {}
        virtual double payoff(double price_at_maturity) const = 0;
    };

    struct european_call: european {
        inline european_call(double strike, datetime const& maturity): european{strike, maturity} {}
        inline european_call(double strike, year_month_day const& maturity): european{strike, maturity} {}
        inline double payoff(double price_at_maturity) const override {
            return std::max(price_at_maturity-K,0.0);
        }
    };

    struct european_put: european {
        inline european_put(double strike, datetime const& maturity): european{strike, maturity} {}
        inline european_put(double strike, year_month_day const& maturity): european{strike, maturity} {}
        inline double payoff(double price_at_maturity) const override {
            return std::max(K-price_at_maturity,0.0);
        }
    };

    struct american: instrument {
        inline american(double strike, datetime const& maturity): instrument{strike,maturity} {} // K{strike}, maturity{maturity} {}
        virtual double payoff(double price) const = 0;
    };

    struct american_call: american {
        american_call(double strike, datetime const& maturity): american{strike, maturity} {}
        american_call(double strike, year_month_day const& maturity): american{strike, maturity} {}
        inline double payoff(double price) const override {
            return std::max(price - K, 0.0);
        }
    };

    struct american_put: american {
        american_put(double strike, datetime const& maturity): american{strike, maturity} {}
        american_put(double strike, year_month_day const& maturity): american{strike, maturity} {}
        inline double payoff(double price) const override {
            return std::max(K - price, 0.0);
        }
    };

}

#endif //BSM_INSTRUMENTS_H
