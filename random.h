#ifndef BLACKSCHOLES_RANDOM_H
#define BLACKSCHOLES_RANDOM_H
#include <random>
#include <ranges>
#include <algorithm>
#include <vector>

/**
 * Random Normal Generator
 * @tparam T
 */
template<typename T>
class random_normal {
private:
    std::random_device device{};
    std::mt19937 generator{device()};
    std::normal_distribution<T> normal;
public:
    inline random_normal() = default;
    inline T operator()() {
        return normal(generator);
    }
    std::vector<T> operator()(int n) {
        std::vector<T> z(n);
        std::ranges::generate(z, [this] () { return this->operator()(); } );
        return z;
    }
};

#endif //BLACKSCHOLES_RANDOM_H
