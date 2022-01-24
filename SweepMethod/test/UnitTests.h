#pragma once

#include <functional>
#include <iostream>
#include <vector>

#include "TestRunner.h"
#include "Profiler.h"

#include "main/interfaces/Instrumental.h"
#include "main/interfaces/serial/SerialInstrumental.h"
#include "main/interfaces/serial/SerialSweepMethod.h"
#include "main/interfaces/parallel/ParallelInstrumental.h"
#include "main/interfaces/parallel/ParallelSweepMethod.h"


class UnitTests {
private:
    size_t N, node;
    vec v, u;

    double h;
    vec x;
    vec A, C, B;

    vec Phi;
    pairs kappa, mu, gamma;

    size_t threadNum, blockSize, classicSize;
    matr pA;
    vec pB;

    void prepareSerialDataForTest(const SerialSweepMethod& sweepMethod) {
        std::tie(v, u,
                 N, node, h, x,
                 A, C, B,
                 Phi, kappa, mu, gamma) = sweepMethod.getFields();
    }

    void setSerialFields(SerialSweepMethod& sweepMethod) {
        sweepMethod.setAllFields(v, u,
                                 N, node, h, x,
                                 A, C, B,
                                 Phi, kappa, mu, gamma);
    }

    void prepareParallelDataForTest(const ParallelSweepMethod& sweepMethod) {
        std::tie(N,
                 threadNum, blockSize, classicSize,
                 pA, pB) = sweepMethod.getAllFields();
    }

    void setParallelFields(ParallelSweepMethod& sweepMethod) {
        sweepMethod.setAllFields(N, threadNum, blockSize, classicSize,pA, pB);
    }

public:
    static void testExecutionTime() {
        std::cout << "I) Execution time for parallel:\n";
        {
            LOG_DURATION("3000 N, 6 threadNum");

            ParallelInstrumental pi(3000, 6);
            matr a = pi.createThirdDiagMatrI();
        }

        {
            LOG_DURATION("3000 N, 1 threadNum");

            ParallelInstrumental pi(3000, 1);
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

        {
            ParallelInstrumental pi(10000, 10);
            ASSERT(pi.checkData());
        }
    }

    void testModelTask() {
        std::cout << "III) Module 7. Model task for estimating computational error:\n";

        {
            SerialSweepMethod ssm(5);
            this->prepareSerialDataForTest(ssm);

            ASSERT(N > 1);

            for (size_t i = 0; i < node; i++) {
                u[i] = 10 + 90 * x[i] * x[i];
            }

            double total = 12. / (h * h);
            A.assign(N, total);
            C.assign(N, 2. * total + 5.);
            B.assign(N, total);

            for (size_t i = 0; i < node - 2; i++) {
                Phi[i] = 450. * x[i + 1] * x[i + 1] - 2110.;
            }

            mu = std::make_pair(10., 100.);
            kappa = std::make_pair(0., 0.);
            gamma = std::make_pair(1., 1.);

            this->setSerialFields(ssm);

            v = ssm.run();
            ASSERT_FOR_DOUBLES(u, v);

            Phi.assign(node, 1);
            Phi[0] = mu.first; Phi[node - 1] = mu.second;
            for (size_t i = 1; i < node - 1; i++) {
                Phi[i] = 2110. - 450. * x[i] * x[i];
            }

            this->setSerialFields(ssm);

            vec res = Instrumental::calcMatrVecMult(ssm.createMatr(), v);
            printf("The scheme (SLAU) is solved with a discrepancy ||R|| = %f\n", ssm.calcR(res, Phi));
            printf("Estimation of the scheme error Z = %f\n", ssm.calcZ());
        }
    }

    void testParallelAlgorithm() {
        std::cout << "IV) Checking the correctness of the parallel algorithm:\n";

        {
            ParallelSweepMethod psm(16, 4);
            this->prepareParallelDataForTest(psm);

            ASSERT(N > 1);
        }
    }

    void execute() {
        TestRunner testRunner;
        str line = "-------------------------------";

        std::vector<std::function<void()>> tests = {
            []() { UnitTests::testExecutionTime(); },
            []() { UnitTests::testEnteredData(); },
            [this]() { this->testModelTask(); },
            [this]() { this->testParallelAlgorithm(); }
        };

        for (const auto& test : tests) {
            std::cout << line << std::endl;
            RUN_TEST(testRunner, test);
        }

        std::cout << line << std::endl;
    }
};