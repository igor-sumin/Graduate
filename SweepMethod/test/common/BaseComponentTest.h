#pragma once

#include <iostream>
#include <vector>
#include <tuple>
#include <random>
#include <functional>
#include <algorithm>
#include <numeric>

#include "main/interfaces/Instrumental.h"

class BaseComponentTest {
protected:
    size_t N, node;
    vec v, u;

public:
    static void execute(const std::vector<std::function<void()>>& tests) {
        TestRunner testRunner;
        str line = "-------------------------------";

        for (const auto& test : tests) {
            std::cout << line << std::endl;
            RUN_TEST(testRunner, test);
        }

        std::cout << line << std::endl;
    }
};