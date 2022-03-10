#pragma once

class MeasureTimeComponentTest : public BaseComponentTest {
public:
    MeasureTimeComponentTest() : BaseComponentTest() {}

//    void testExecutionTimeBetweenAlgorithms() {
//
//        Phi.assign(node, 1);
//        Phi[0] = mu.first; Phi[node - 1] = mu.second;
//        for (size_t i = 1; i < node - 1; i++) {
//            Phi[i] = 2110. - 450. * x[i] * x[i];
//        }
//    }

    static void testCompareAlgorithmsByTask7(int n, int tN) {
        vec ser, par;

        {
            LOG_DURATION("Serial for " + std::to_string(n))

            SerialSweepMethod ssm(n, TASK::TASK_7);
            vec res = ssm.run();
            Instrumental::printVec(res, "serial answer");

            ser = res;
        }

        {
            LOG_DURATION("Parallel for " + std::to_string(n))

            ParallelSweepMethod psm(n, tN, TASK::TASK_7);
            vec res = psm.run();
            Instrumental::printVec(res, "parallel answer");

            par = res;
        }
    }

    void execute() {
        std::vector<std::function<void()>> tests = {
            []() { MeasureTimeComponentTest::testCompareAlgorithmsByTask7(16, 4); }
        };

        BaseComponentTest::execute(tests, "Measure Time Component Test");
    }
};