#  BSM

This is a personal project created to learn c++ 20 by implementing option pricing and greeks calculation under the Black Scholes Merton model. You may find bugs, convoluted design and some horrible copy-paste. This code is an experimentation, do not expect any good coding practices or idiomatic c++. I've thrown away 2 other iterations of this same codebase and it might happen again. 

See the CMakeLists.txt for the required dependencies. I run Win11 and used GCC 11.2 on WSL 2.0, and everything was easy to install.

Designed and implemented using TDD/Catch2 while also performing micro benchmark using Google Benchmark. The code has not been optimized yet. 
Focus has been on hiding implementation layer complexities while providing expressiveness to library users. It uses the nice Autodiff library (https://autodiff.github.io).

### TODO

* Implement methods that are actually more interesting, not the closed formulas for european options.
* Code refactoring, optimization, etc.
* More benchmarks

### NOTES

* Started implementing the QD+ approximation method for American Options, but there is a lot of work to be done. So far it works for american puts, and I've written two tests based on the data of the original paper. There is another test which shows it coincides with the binomial tree result too.
* Also implemented a basic binomial tree method for easy of comparison (using Cox-Ross-Rubinstein method). But I need to improve it to show the exercise boundary and consolidate the interface used for american option pricing.

## Example 

    auto K = 100.0;
    auto S = 100.0;
    auto sigma = 0.20;
    auto t = system_clock::now();
    auto r = 0.01;
    auto q = 0.05;
    mkt_params mktParams{S, sigma, t, r, q};
    european_call europeanCall{K, t + 0.5_years};
    analytical_solver<autodiff_off> solve1{mktParams};
    analytical_solver<autodiff_dual> solve2{mktParams};
    analytical_solver<autodiff_var> solve3{mktParams};

    auto callPricing = solve1(europeanCall);
    auto callPricing_dual = solve2(europeanCall);
    auto callPricing_var = solve3(europeanCall);

    //Check Prices against the analytical formula
    CHECK(callPricing_dual->price()==Approx(callPricing->price()));
    CHECK(callPricing_var->price()==Approx(callPricing->price()));

    //Check Prices
    CHECK(callPricing_dual->price()==Approx(callPricing_var->price()));

    //Check Delta
    CHECK(callPricing_dual->delta()==Approx(callPricing_var->delta()));

    //Check Gamma
    CHECK(callPricing_dual->gamma()==Approx(callPricing_var->gamma()));