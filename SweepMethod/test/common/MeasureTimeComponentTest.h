#pragma once

class MeasureTimeComponentTest : public BaseComponentTest {
public:
    MeasureTimeComponentTest() : BaseComponentTest() {}

    static void testCompareAlgorithmsByTask7(int n, int tN) {
        vec ser, par;

        {
            LOG_DURATION("Serial for " + std::to_string(n))

            SerialSweepMethod ssm(n, Instrumental::TASK::TASK_7);
            vec res = ssm.run();
            Instrumental::printVec(res, "serial answer");

            ser = res;
        }

        {
            LOG_DURATION("Parallel for " + std::to_string(n))

            ParallelSweepMethod psm(n, tN, Instrumental::TASK::TASK_7);
            vec res = psm.run();
            Instrumental::printVec(res, "parallel answer");

            par = res;
        }

        printf("compare = %d", Instrumental::compareVec(ser, par));
    }

    void execute() {
        std::vector<std::function<void()>> tests = {
            []() { MeasureTimeComponentTest::testCompareAlgorithmsByTask7(9, 3); }
        };

        BaseComponentTest::execute(tests, "Measure Time Component Test");
    }
};