#ifndef BLACKSCHOLES_BSM_H
#define BLACKSCHOLES_BSM_H

#include <autodiff/reverse/var.hpp>

using namespace autodiff;

enum class CP { call, put };
const double one_div_root_two = 1.0/sqrt(2.0);
const double one_div_root_two_pi = 1.0/sqrt(2.0*std::numbers::pi);

/**
 * Standard Normal CDF
 * @param x
 * @return
 */
var phi(var const& x);

/**
 * Calculates the premium of a call or put european option
 * @param cp Defines if function will price a call or a put
 * @param K Option Strike
 * @param S Underlying Spot Price
 * @param variance Variance (square of the sigma)
 * @param tau Time to maturity
 * @param r Discounting rate
 * @param q Dividend yield (aka Convenience yield)
 * @return Premium of Call or put
 */
var european(CP cp, double K, var const& S, var const& variance, var const& tau, var const& r, var const& q);

/**
 * Calculates price of a Forward, ie, the strike for a 0-valued forward.
 * @param S Instrument spot Price
 * @param tau Time to maturity (1 = 1 year)
 * @param r Discounting rate
 * @param q Dividend yield (aka Convenience yield)
 * @return Price of a Forward (not the value)
 */
var forward(var const& S, var const& tau, var const& r, var const& q);

/**
 * Calculates the value of a Forward with strike K.
 * @param K Forward Strike
 * @param S Instrument spot Price
 * @param tau Time to maturity (1 = 1 year)
 * @param r Discounting rate
 * @param q Dividend yield (aka Convenience yield)
 * @return Value of the Forward
 */
var forward(double K, var const& S, var const& tau, var const& r, var const& q);

/**
 * Calculates the implied dividend using Newton-Raphson method.
 * The derivative was calculated using reverse-mode autodiff.
 * @param K
 * @param F
 * @param S
 * @param tau
 * @param r
 * @param q Starting value for Implied Dividend search
 * @return Implied Dividend yield
 */
var implied_dividend(double K, var const& F, var const& S, var const& tau, var const& r, var q = 0.0);

/**
 * Calculates the implied volatility using Newton-Raphson method.
 * The derivative was calculated using reverse-mode autodiff.
 * @param cp
 * @param K
 * @param P
 * @param S
 * @param tau
 * @param r
 * @param q
 * @param sigma Starting value for Implied Vol search
 * @return Implied Vol
 */
var implied_volatility(CP cp, double K, var const& P, var const& S, var const& tau, var const& r, var const& q, var sigma = 0.10);

#endif //BLACKSCHOLES_BSM_H
