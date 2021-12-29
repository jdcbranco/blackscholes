#include <autodiff/reverse/var.hpp>
#include <autodiff/reverse/var/eigen.hpp>
#include "options.h"

using namespace autodiff;

const double one_div_root_two = 1.0/sqrt(2.0);
const double one_div_root_two_pi = 1.0/sqrt(2.0*std::numbers::pi);

namespace options {

    template<typename T>
    T phi(T const &x) {
        return 0.5 * (1.0 + erf(x * one_div_root_two));
    }

    template<>
    var phi(var const &x) {
        return 0.5 * (1.0 + erf(x * one_div_root_two));
    }

    template<>
    std::shared_ptr<detail::Expr<double>> phi(std::shared_ptr<detail::Expr<double>> const &x) {
        return 0.5 * (1.0 + erf(x * one_div_root_two));
    }

}