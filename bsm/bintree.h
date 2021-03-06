#ifndef BSM_BINTREE_H
#define BSM_BINTREE_H

#include <cassert>
#include <vector>
#include <iostream>
#include <iomanip>

namespace bsm {
    template<typename T>
    class bintree {

    public:
        explicit bintree(int steps) : steps_(steps), lattice(steps) {
            for (int i = 0; i < steps; ++i) {
                lattice[i].resize(i + 1);
            }
        }

        bintree(bintree const &copy) = default;

        T operator()(int t, int i) const {
            assert(("Invalid index t", t >= 0 and t < steps_));
            assert(("Invalid index i", i >= 0 and i <= t));
            return lattice[t][i];
        }

        auto operator()(int t) {
            assert(("Invalid index t", t >= 0 and t < steps_));
            return make_tuple(lattice[t].begin(), lattice[t].end());
        }

        void set(int t, int i, T value) {
            assert(("Invalid index t", t >= 0 and t < steps_));
            assert(("Invalid index i", i >= 0 and i <= t));
            lattice[t][i] = value;
        }

        int size() const { return steps_; }

        T root() const { return lattice[0][0]; }

    private:
        int steps_;
        std::vector<std::vector<T>> lattice;
    };

    template<typename T>
    std::ostream &operator<<(std::ostream &out, bintree<T> const &tree) {
        int size = tree.size();
        int tabs = 1 + (size + 1) / 2;

        for (int t = 0; t < size; ++t) {
            for (int k = 0; k < tabs; ++k) {
                out << std::setw(4) << "\t";
            }
            tabs -= 1;
            for (int i = 0; i <= t; ++i) {
                T value = tree(t, i);
                out << value << std::setw(4) << "\t";
            }
            out << "\n";
        }
        return out;
    }
}
#endif //BSM_BINTREE_H
