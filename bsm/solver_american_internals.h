#ifndef BSM_SOLVER_AMERICAN_INTERNALS_H
#define BSM_SOLVER_AMERICAN_INTERNALS_H

#include "solver.h"

#include <cassert>

namespace bsm {
    namespace internals {

        template<typename T>
        inline bool never_optimal_exercise(pricing<T> const& p, bool call) {
            return call ? (p.q<=0 and (p.q <= p.r)) : (p.r <= 0 and (p.r <= p.q));
        }

        template<typename T>
        inline bool never_optimal_exercise(pricing<T> const& p, instrument_type const& type) {
            assert(("Haven't implemented optimal exercise for non-options", type==instrument_type::call or type==instrument_type::put));
            return call ? (p.q<=0 and (p.q <= p.r)) : (p.r <= 0 and (p.r <= p.q));
        }

        template<typename T>
        inline T exercise_boundary_at_maturity(pricing<T> const& p, instrument_type const& type) {
            switch (type) {
                case instrument_type::call:
                        if(p.r<=p.q) {
                        return p.K;
                        } else {
                        return p.K * p.r/p.q;
                        }
                case instrument_type::put:
                        if(p.r>=p.q) {
                        return p.K;
                        } else {
                        return p.K * p.r/p.q;
                        }
                default:
                    return p.K;
            }
        }

    }
}

#endif //BSM_SOLVER_AMERICAN_INTERNALS_H
