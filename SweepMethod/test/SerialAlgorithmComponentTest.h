#pragma once

#include <functional>
#include <iostream>
#include <vector>
#include <test/common/BaseComponentTest.h>


class SerialAlgorithmComponentTest final : public SerialSweepMethod, public BaseComponentTest {
public:
    SerialAlgorithmComponentTest() : SerialSweepMethod() {}

    /*
     * Module 7. Model task for estimating computational error (Serial realization)
     */
    static void testTask7(size_t n) {
        SerialSweepMethod ssm(n, Task::TASK::TASK_7);
        vec ar = ssm.run();
        vec er = ssm.createResByTask7();

        Instrumental::printVec(ar, "ar");
        Instrumental::printVec(er, "er");

        Instrumental::compareVec(ar, er);
    }

    static void execute() {
        std::vector<std::function<void()>> tests = {
                []() { SerialAlgorithmComponentTest::testTask7(10); }
        };

        BaseComponentTest::execute(tests, "Serial Component Test");
    }
};