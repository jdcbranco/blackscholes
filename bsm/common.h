#ifndef BSM_COMMON_H
#define BSM_COMMON_H

#include <concepts>
#include <utility>
#include "chrono.h"

namespace bsm {
    using namespace bsm::chrono;

    template<typename T = double> //requires std::floating_point<T> //, typename CLOCK = system_clock
    struct mkt_params {
        T S;
        T sigma;
        datetime t;
        T r;
        T q;

        mkt_params(T const& S, T const& sigma, datetime const& t, T const& r, T const& q):
                S{S}, sigma{sigma}, t{t}, r{r}, q{q}
        {
        }

        mkt_params(T const& S, T const& sigma, std::chrono::year_month_day const& t, T const& r, T const& q):
                S{S}, sigma{sigma}, t{t}, r{r}, q{q}
        {
        }

        mkt_params(mkt_params<long double> const& mp): mkt_params{mp.S, mp.sigma, mp.t, mp.r, mp.q} {}

        template<std::size_t Index>
        std::tuple_element_t<Index, mkt_params<T>>& get() &
        {
            if constexpr (Index == 0) return S;
            if constexpr (Index == 1) return sigma;
            if constexpr (Index == 2) return t;
            if constexpr (Index == 3) return r;
            if constexpr (Index == 4) return q;
        }

        template<std::size_t Index>
        std::tuple_element_t<Index, mkt_params<T>> const& get() const&
        {
            if constexpr (Index == 0) return S;
            if constexpr (Index == 1) return sigma;
            if constexpr (Index == 2) return t;
            if constexpr (Index == 3) return r;
            if constexpr (Index == 4) return q;
        }

        template<std::size_t Index>
        std::tuple_element_t<Index, mkt_params<T>>& get() &&
        {
            if constexpr (Index == 0) return std::move(S);
            if constexpr (Index == 1) return std::move(sigma);
            if constexpr (Index == 2) return std::move(t);
            if constexpr (Index == 3) return std::move(r);
            if constexpr (Index == 4) return std::move(q);
        }

        template<std::size_t Index>
        std::tuple_element_t<Index, mkt_params<T>> const& get() const&&
        {
            if constexpr (Index == 0) return std::move(S);
            if constexpr (Index == 1) return std::move(sigma);
            if constexpr (Index == 2) return std::move(t);
            if constexpr (Index == 3) return std::move(r);
            if constexpr (Index == 4) return std::move(q);
        }
    };
}

namespace std {
    template<typename T>
    struct tuple_size<bsm::mkt_params<T>> : std::integral_constant<size_t, 5> { };

    template<size_t Index, typename T>
    struct tuple_element<Index, bsm::mkt_params<T>>
    : tuple_element<Index, tuple<T, T, bsm::chrono::datetime, T, T>>
    {
    };
}


#endif //BSM_COMMON_H
