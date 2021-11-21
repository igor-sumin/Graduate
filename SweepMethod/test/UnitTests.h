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

public:
	static void testEnteredData() {
		std::cout << "Checking the correctness of the entered data\n";

		{
			ParallelInstrumental pi(8, 4, -1, -1);
			ASSERT(!pi.checkData());
		}
		{
			ParallelInstrumental pi(80, 10, -1, -1);
			ASSERT(!pi.checkData());
		}
		{
			ParallelInstrumental pi(100, 10, -1, -1);
			ASSERT(!pi.checkData());
		}
	}

	/*
	void test1() {
		this->prepareDataForTest();

		// Sequential implementation of the sweep method
		{
			LOG_DURATION("SerialSweepMethod");

			pairs kappa = std::make_pair(0.5, 0.5);
			pairs mu = std::make_pair(0, N - 1);
			pairs gamma = std::make_pair(1., 1.);

			vec a(N - 2, 1.);
			vec c(N - 2, -3.);
			vec b(N - 2, 2.);
			vec f(N - 2);
			std::iota(f.begin(), f.end(), 1.);
			std::for_each(f.begin(), f.end(), [](double& n) { n *= -1; });

			SerialSweepMethod ssm(a, c, b, f, kappa, mu, gamma);
			vec x = ssm.run();

			Instrumental::printVec(x, "\nSerial x");

			std::cout << "The scheme (SLAU) is solved with a discrepancy ||R|| = "
				<< I.getR(I.calcMatrVecMult(A, x), B) << std::endl;
		}

		std::cout << std::endl;

		// Parallel implementation of the run
		{
			LOG_DURATION("ParallelSweepMethod");

			ParallelSweepMethod psm(A, B);
			vec x = psm.run();

			Instrumental::printVec(x, "\nParallel x");

			std::cout << "The scheme (SLAU) is solved with a discrepancy ||R|| = "
				<< I.getR(I.calcMatrVecMult(A, x), B) << std::endl;
		}
	}

	void test2() {
		this->prepareDataForTest();

		// Sequential implementation of the run (with an exact solution)
		{
			LOG_DURATION("SerialSweepMethod");

			vec xi(N);
			double h = 1 / N;
			for (int i = 0; i < N; i++) {
				xi[i] = i * h;
			}

			/*
			std::cout << "The scheme (SLAU) is solved with a discrepancy ||R|| = " << I.getR(calcMatrVecMult(A, x), B) << std::endl;
			std::cout << "Estimation of the scheme error Z = " << I.getZ(u, x) << std::endl;
		}
	}
	*/

	void testModelTask() {
		std::cout << "Module 7. Model task for estimating computational error\n";

		{
            SerialSweepMethod ssm(5);
            std::tie(v, u, N, node, h, x, A, C, B, Phi, kappa, mu, gamma) = ssm.getAllFields();
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

            ssm.setAllFields(v, u, N, node, h, x, A, C, B, Phi, kappa, mu, gamma);

            v = ssm.run();
            ASSERT_FOR_DOUBLES(u, v);

            Phi.assign(node, 1);
            Phi[0] = mu.first; Phi[node - 1] = mu.second;
            for (size_t i = 1; i < node - 1; i++) {
                Phi[i] = 2110. - 450. * x[i] * x[i];
            }

            ssm.setAllFields(v, u, N, node, h, x, A, C, B, Phi, kappa, mu, gamma);

            vec res = Instrumental::calcMatrVecMult(ssm.createMatr(), v);
            printf("The scheme (SLAU) is solved with a discrepancy ||R|| = %f\n", ssm.calcR(res, Phi));
            printf("Estimation of the scheme error Z = %f\n", ssm.calcZ());
		}
	}

	void execute() {
		TestRunner testRunner;
		str line = "-------------------------------";

		std::vector<std::function<void()>> tests = {
			[this]() { this->testModelTask(); }
		};

		for (const auto& test : tests) {
			std::cout << line << std::endl;

			RUN_TEST(testRunner, test);

			std::cout << line << std::endl;
		}
	}
};