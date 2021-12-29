#  BSM Calculations in C++ using Autodiff

This is a personal project created to learn c++ 17/20 while trying the new Autodiff c++ library (https://autodiff.github.io) in order to compute option pricing and greeks calculation. 
The expressiveness of this new library made me quite impressed. However, the reverse mode so far consumes too much memory to be worthwhile in Monte Carlo method or even Binomial trees, and thus so far I only managed to use it for simple European-style options (call, puts) and Forwards with straightforward closed-form analytical expressions. 

See the CMakeLists.txt for the required dependencies. I run Win11 and used GCC 11.2 on WSL 2.0, and everything was easy to install.

There are some sample unit tests (using Catch2) and sample benchmark (using Google Benchmark). The code has not been optimized yet, and I only tried the Reverse-mode AD so far. 

### TODO
* Try the Forward-mode and compare performance with Reverse-mode
* Implement more methods American Option Pricing methods
* Code refactoring, optimization, etc.

## Example 

    //Example of European Call pricing and greek calculation using autodiff
    double const K = 15000.0;

    options::params<var> main_params{15000.0, 0.21, 0.5, 0.0, 0.0};
    auto&& [S, sigma, tau, r, q] = main_params;

    options::european_call<var> europeanCall{K};
    analytical::analyticalmethod<var> bsm_reverse_mode{main_params};

    auto call = bsm_reverse_mode.solve(europeanCall);
    auto [delta, gamma, vega, theta, rho, psi] = all_greeks(call); // greeks;

    std::cout << "call = " << call << std::endl;
    std::cout << "delta = " << delta << std::endl;
    std::cout << "gamma = " << gamma << std::endl;
    std::cout << "vega = " << vega << std::endl;
    std::cout << "theta = " << theta << std::endl;
    std::cout << "rho = " << rho << std::endl;
    std::cout << "psi = " << psi << std::endl;

    //Example of Implied Dividend calculation using Newton-Raphson method and autodiff
    options::forward<var> fwdInstrument{K};
    auto fwd = bsm_reverse_mode.solve(fwdInstrument);
    auto impl_div = bsm_reverse_mode.imply_dividends(fwdInstrument, (var)fwd);

    std::cout << "fwd = " << fwd << std::endl;
    std::cout << std::boolalpha;
    std::cout << "impl div = " << impl_div << " which is " << (abs(impl_div - q) <= 1.0e-4) << std::endl;

    //Example of Implied Volatility calculation using Newton-Raphson method and autodiff
    auto iv = bsm_reverse_mode.imply_volatility(europeanCall, val((var)call));

    std::cout << "impl vol = " << iv << " which is " << (abs(iv - sigma) <= 1.0e-9) << std::endl;
