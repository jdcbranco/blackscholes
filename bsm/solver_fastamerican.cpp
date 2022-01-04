#include "solver.h"

#include <vector>
#include <algorithm>

namespace bsm {

    struct fastamerican_method: pricing<double> {
        const int l;
        const int m;
        const int n;
        const double tau_max;

        fastamerican_method(american const& instrument, mkt_params<double> mp, int l, int m, int n):
        pricing{instrument,mp},
        l{l}, m{m}, n{n}, tau_max{tau}
        {

        }

        std::vector<double> get_chebyshev_nodes(int _n) {
            std::vector<double> z(_n+1);
            std::generate(z.begin(), z.end(), [i = 0,_n] () mutable { return -cos(i++*std::numbers::pi/_n); });
            return z;
        }


    };


}