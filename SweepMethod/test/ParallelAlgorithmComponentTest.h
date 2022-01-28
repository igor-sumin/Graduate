#pragma once

#include <main/interfaces/parallel/ParallelSweepMethod.h>
#include <test/common/Profiler.h>
#include <test/common/TestRunner.h>
#include <test/common/BaseComponentTest.h>

class ParallelAlgorithmComponentTest final : public ParallelSweepMethod, public BaseComponentTest {
private:
    void prepareParallelDataForTest(const ParallelSweepMethod& sweepMethod) {
        std::tie(N,
                 threadNum, blockSize, interSize,
                 A, b, y) = sweepMethod.getAllFields();
    }

    void setParallelFields(ParallelSweepMethod& sweepMethod) {
        sweepMethod.setAllFields(N, threadNum, blockSize, interSize, A, b, y);
    }

public:
    ParallelAlgorithmComponentTest() : ParallelSweepMethod() {}

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
                {1.000,    0.000,    0.000,    0.000,    0.000,    0.000},
                {300.000, -605.000,  300.000,  0.000,    0.000,    0.000},
                {0.000,    300.000, -605.000,  300.000,  0.000,    0.000},
                {0.000,    0.000,    300.000, -605.000,  300.000,  0.000},
                {0.000,    0.000,    0.000,    300.000, -605.000,  300.000},
                {0.000,    0.000,    0.000,    0.000,    0.000,    1.000}
            };
            b = {10.0, 2092.0, 2038.0, 1948.0, 1822.0, 100.0};
            u = {10, 13.6, 24.4, 42.4, 67.6, 100};
            this->setParallelFields(psm);

            // execution
            psm.transformation();
            this->prepareParallelDataForTest(psm);

            // check
            y.resize(A.size());
            for (int i = 0; i < A.size(); i++) {
                y[i] = b[i] / A[i][i];
            }

            printf("The scheme (SLAU) is solved with a discrepancy ||R|| = %f\n", si.calcR(y, u));
        }

        print();

        {
            LOG_DURATION("time (8, 2)");
            ParallelSweepMethod psm(8, 2);

            psm.transformation();
            this->prepareParallelDataForTest(psm);

            Instrumental::printMatr(A, "A (8, 2)");
            Instrumental::printVec(b, "b (8, 2)");
        }

        print();

        {
            LOG_DURATION("time (12, 3)");
            ParallelSweepMethod psm(12, 3);

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

    std::pair<matr, vec> testCollectInterferElemPreprocessing(int n, int tN) {
        ParallelSweepMethod psm(n, tN);

        psm.transformation();
        this->prepareParallelDataForTest(psm);

        {
            LOG_DURATION("serial");

            matr R(interSize, vec(interSize, 0.));
            vec partB(interSize);

            size_t bS1 = blockSize - 1;
            R[0][0] = A[bS1][bS1];
            R[0][1] = A[bS1][blockSize];
            R[1][0] = A[blockSize][bS1];
            R[2][1] = A[2 * blockSize - 1][bS1];

            size_t bSN = N - blockSize;
            R[interSize - 1][interSize - 1] = A[bSN][bSN];
            R[interSize - 1][interSize - 2] = A[bSN][bSN - 1];
            R[interSize - 2][interSize - 1] = A[bSN - 1][bSN];
            R[interSize - 3][interSize - 2] = A[N - blockSize * 2][bSN];

            for (size_t i = blockSize, k = 1; i < N; i += blockSize, k += 2) {
                partB[k - 1] = b[i - 1];
                partB[k] = b[i];
            }

            // printMatr(R, "R11");
        }

        print();

        {
            LOG_DURATION("parallel");

            matr R(interSize, vec(interSize, 0.));
            vec partB(interSize);

            size_t iter;
            size_t s;

            #pragma omp parallel private(s, iter) shared(R, partB) num_threads(threadNum - 1) default(none)
            {
                iter = omp_get_thread_num();

                #pragma omp taskgroup
                {
                    // (extreme parts 1)
                    #pragma omp task shared(R) default(none)
                    this->preULR(R);

                    // (extreme parts 2)
                    #pragma omp task shared(R) default(none)
                    this->preLRR(R);

                    // filling intersecting elements from vector b to partB
                    #pragma omp taskloop private(s) firstprivate(iter) shared(partB) default(none)
                    for (s = iter; s < iter + 1; s++) {
                        partB[2 * s] = b[(s + 1) * blockSize - 1];
                        partB[2 * s + 1] = b[(s + 1) * blockSize];
                    }
                }
            }

            // printMatr(R, "R21");

            return std::make_pair(R, partB);
        }
    }

    std::pair<matr, vec> testCollectInterferElemPostprocessing(int n, int tN) {
        ParallelSweepMethod psm(n, tN);
        matr R; vec partB;

        psm.transformation();

        {
            LOG_DURATION("serial");

            std::tie(R, partB) = this->testCollectInterferElemPreprocessing(n, tN);
            this->prepareParallelDataForTest(psm);

            // 2. post-processing (internal part)
            for (size_t i = blockSize, k = 1; i < N - blockSize; i += blockSize, k += 2) {
                for (size_t j = blockSize, l = 1; j < N - blockSize; j += blockSize, l += 2) {

                    // a1 ----- a2
                    // |         |
                    // |         |
                    // a3 ----- a4
                    double a1 = A[i][j];
                    double a2 = A[i][j + blockSize - 1];
                    double a3 = A[i + blockSize - 1][j];
                    double a4 = A[i + blockSize - 1][j + blockSize - 1];

                    if (a1 != 0 && a4 != 0) {
                        R[k][l] = a1;
                        R[k + 1][l + 1] = a4;
                    } else if (a1 != 0) {
                        R[k][l - 1] = a1;
                        R[k + 1][l] = a3;
                    } else if (a4 != 0) {
                        R[k][l + 1] = a2;
                        R[k + 1][l + 2] = a4;
                    }
                }
            }

            // printMatr(R, "R12");
        }

        print();

        {
            LOG_DURATION("parallel");

            std::tie(R, partB) = this->testCollectInterferElemPreprocessing(n, tN);
            this->prepareParallelDataForTest(psm);

            size_t k = 1, l = 1;
            size_t iter, jter;
            size_t i, j;
            double a1, a2, a3, a4;

            #pragma omp parallel private(iter, jter, i, j, a1, a2, a3, a4) firstprivate(k, l) shared(A, R, blockSize) default(none)
            {
                iter = (omp_get_thread_num() + 1) * blockSize;

                for (i = iter; i < N - iter; i += iter) {
                    for (j = iter; j < N - iter; j += iter) {
                        // a1 ----- a2
                        // |         |
                        // |         |
                        // a3 ----- a4
                        a1 = A[i][j];
                        a2 = A[i][j + blockSize - 1];
                        a3 = A[i + blockSize - 1][j];
                        a4 = A[i + blockSize - 1][j + blockSize - 1];

                        if (a1 != 0 && a4 != 0) {
                            R[k][l] = a1;
                            R[k + 1][l + 1] = a4;
                        } else if (a1 != 0) {
                            R[k][l - 1] = a1;
                            R[k + 1][l] = a3;
                        } else if (a4 != 0) {
                            R[k][l + 1] = a2;
                            R[k + 1][l + 2] = a4;
                        }

                        l += 2;
                    }

                    l = 1;
                    k += 2;
                }
            }

            // printMatr(R, "R22");

            return std::make_pair(R, partB);
        }
    }

    void testOrderingCoefficient(int n, int tN) {
        ParallelSweepMethod psm(n, tN);
        matr R; vec partB;

        psm.transformation();

        {
            LOG_DURATION("serial");

            std::tie(R, partB) = this->testCollectInterferElemPostprocessing(n, tN);
            this->prepareParallelDataForTest(psm);

            for (int i = 0; i < interSize; i += 2) {
                std::swap(R[i][i], R[i][i + 1]);
                std::swap(R[i + 1][i], R[i + 1][i + 1]);
                std::swap(partB[i], partB[i + 1]);
            }

            // printMatr(R, "R3");
            // printVec(partB, "partB3");
        }

        {
            LOG_DURATION("parallel");

            std::tie(R, partB) = this->testCollectInterferElemPostprocessing(n, tN);
            this->prepareParallelDataForTest(psm);

            #pragma omp parallel for shared(R, partB) default(none)
            for (int i = 0; i < interSize; i += 2) {
                std::swap(R[i][i], R[i][i + 1]);
                std::swap(R[i + 1][i], R[i + 1][i + 1]);
                std::swap(partB[i], partB[i + 1]);
            }

            // printMatr(R, "R4");
            // printVec(partB, "partB4");
        }
    }

    void testSlide11() {
        ParallelSweepMethod psm(12, 4);

        psm.transformation();
        psm.testing();

        // std::tie(R, y1) = psm.collectInterferElem();

        this->prepareParallelDataForTest(psm);

//        Instrumental::printMatr(A, "A (12, 3)");
//        Instrumental::printVec(b, "b (12, 3)");
    }

    void execute() {
        std::vector<std::function<void()>> tests = {
//            []() { ParallelAlgorithmComponentTest::testExecutionTime(); },
//            []() { ParallelAlgorithmComponentTest::testEnteredData(); },
//            [this]() { this->testSlide3(); }
//            [this]() { this->testSlide11(); }
//            [this]() { this->testCollectInterferElemPreprocessing(12, 4); },
//            [this]() { this->testCollectInterferElemPostprocessing(12, 3); }
            [this]() { this->testOrderingCoefficient(12, 4); }
        };

        BaseComponentTest::execute(tests);
    }
};