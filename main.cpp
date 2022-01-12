#include "bsm.h"

#include <iostream>

#include <taskflow/taskflow.hpp>

int main() {
    //See tests for examples

    tf::Executor executor;
    tf::Taskflow taskflow;

    auto [A, B, C, D] = taskflow.emplace(  // create four tasks
            [] () { std::cout << "TaskA\n"; },
            [] () { std::cout << "TaskB\n"; },
            [] () { std::cout << "TaskC\n"; },
            [] () { std::cout << "TaskD\n"; }
    );

    A.precede(B, C);  // A runs before B and C
    D.succeed(B, C);  // D runs after  B and C

    int* vec;
    int first, last;

    auto init = taskflow.emplace([&](){
        first = 0;
        last  = 1000;
        vec = new int[1000];
    });

    auto pf = taskflow.for_each_index(std::ref(first), std::ref(last), 1,
                                      [&] (int i) {
                                          std::cout << "parallel iteration on index " << vec[i] << '\n';
                                      }
    );

    init.precede(pf);

    executor.run(taskflow).wait();

    return 0;
}