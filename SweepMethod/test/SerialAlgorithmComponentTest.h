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
    static void testTask7() {
        SerialSweepMethod ssm(5, Task::TASK::TASK_7);
        vec ar = ssm.run();
        vec er = {10, 13.6, 24.4, 42.4, 67.6, 100};

        Instrumental::compareVec(ar, er) ? printf("ok") : printf("not ok");
    }

    static void execute() {
        std::vector<std::function<void()>> tests = {
                []() { SerialAlgorithmComponentTest::testTask7(); }
        };

        BaseComponentTest::execute(tests, "Serial Component Test");
    }
};