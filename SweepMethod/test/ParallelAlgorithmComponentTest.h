#pragma once

#include <main/interfaces/parallel/ParallelSweepMethod.h>
#include <test/common/Profiler.h>
#include <test/common/TestRunner.h>
#include <test/common/BaseComponentTest.h>

class ParallelAlgorithmComponentTest final : public BaseComponentTest {
private:
    size_t threadNum, blockSize, classicSize;
    matr A;
    vec b, x, u;

    void prepareParallelDataForTest(const ParallelSweepMethod& sweepMethod) {
        std::tie(N,
                 threadNum, blockSize, classicSize,
                 A, b) = sweepMethod.getAllFields();
    }

    void setParallelFields(ParallelSweepMethod& sweepMethod) {
        sweepMethod.setAllFields(N, threadNum, blockSize, classicSize, A, b);
    }

public:
    static void testExecutionTime() {
        std::cout << "I) Execution time for parallel:\n";
        {
            LOG_DURATION("1000 N, 4 threadNum");

            ParallelInstrumental pi(1000, 4);
            matr a = pi.createThirdDiagMatrI();
        }

        {
            LOG_DURATION("10000 N, 1 threadNum");

            ParallelInstrumental pi(1000, 1);
            matr a = pi.createThirdDiagMatrI();
        }

    }

    static void testEnteredData() {
        std::cout << "II) Checking the correctness of the entered data:\n";

        {
            ParallelInstrumental pi(16, 4);
            ASSERT(pi.checkData());
        }

        {
            ParallelInstrumental pi(1600, 10);
            ASSERT(pi.checkData());
        }

        {
            ParallelInstrumental pi(161, 4);
            ASSERT(!pi.checkData());
        }

        {
            ParallelInstrumental pi(4, 4);
            ASSERT(!pi.checkData());
        }

        {
            ParallelInstrumental pi(16, 3);
            ASSERT(!pi.checkData());
        }
    }

    void testSlide3() {
        std::cout << "III) Checking the correctness of the algorithm on the slide 3:\n";

        {
            // creation
            SerialInstrumental si(5);
            ParallelSweepMethod psm(6, 1);
            this->prepareParallelDataForTest(psm);

            A = {
                    {1.000,   0.000,    0.000,    0.000,    0.000,    0.000},
                    {300.000, -605.000, 300.000,  0.000,    0.000,    0.000},
                    {0.000,   300.000,  -605.000, 300.000,  0.000,    0.000},
                    {0.000,   0.000,    300.000,  -605.000, 300.000,  0.000},
                    {0.000,   0.000,    0.000,    300.000,  -605.000, 300.000},
                    {0.000,   0.000,    0.000,    0.000,    0.000,    1.000}
            };
            b = {10.0, 2092.0, 2038.0, 1948.0, 1822.0, 100.0};
            u = {10, 13.6, 24.4, 42.4, 67.6, 100};
            this->setParallelFields(psm);

            // execution
            psm.transformation();
            this->prepareParallelDataForTest(psm);

            // check
            x.resize(A.size());
            for (int i = 0; i < A.size(); i++) {
                x[i] = b[i] / A[i][i];
            }

            printf("The scheme (SLAU) is solved with a discrepancy ||R|| = %f\n", si.calcR(x, u));
        }

        print();

        {
            LOG_DURATION("time (8, 2)");
            ParallelSweepMethod psm(8, 2);
            this->prepareParallelDataForTest(psm);

            A = psm.createThirdDiagMatrI();
            b = psm.createVecN();

            this->setParallelFields(psm);

            psm.transformation();
            this->prepareParallelDataForTest(psm);

            Instrumental::printMatr(A, "A (8, 2)");
            Instrumental::printVec(b, "b (8, 2)");
        }

        print();

        {
            LOG_DURATION("time (12, 3)");
            ParallelSweepMethod psm(12, 3);
            this->prepareParallelDataForTest(psm);

            A = psm.createThirdDiagMatrI();
            b = psm.createVecN();
            this->setParallelFields(psm);

            psm.transformation();
            this->prepareParallelDataForTest(psm);

            Instrumental::printMatr(A, "A (12, 3)");
            Instrumental::printVec(b, "b (12, 3)");
        }

        print();

        {
            // parallel work faster in this case
            //
            // (15000, 3): 9260 ms (9 s)
            // (15000, 1): 14317 ms (14 s)
        }
    }

    void execute() {
        std::vector<std::function<void()>> tests = {
            []() { ParallelAlgorithmComponentTest::testExecutionTime(); },
            []() { ParallelAlgorithmComponentTest::testEnteredData(); },
            [this]() { this->testSlide3(); }
        };

        BaseComponentTest::execute(tests);
    }
};