#pragma once

#include <functional>

#include "TestRunner.h"
#include "Profiler.h"
#include "Instrumental.h"

#include "SerialInstrumental.h"
#include "SerialSweepMethod.h"

#include "ParallelInstrumental.h"
#include "ParallelSweepMethod.h"


class UnitTests {
private:
	size_t N, node; 
	vec u, v;

	int threadNum, blockSize, classicSize;
	matr A;
	vec B;

	double h;
	double a, c, b;
	vec x;

	Instrumental I;
	ParallelInstrumental P;
	SerialInstrumental S;

	void prepareParallelDataForTest() {

		A = P.createThirdDiagMatrI();
		B = P.createVecN();

		tie(N, threadNum, blockSize, classicSize) = P.getAllFields();
	}

	void prepareSerialDataForTest(size_t n) {
		tie(N, node, u, v) = I.getAllFields();
		tie(x, h, a, c, b) = S.getAllFields();
	}

public:
	void testEnteredData() {
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

			std::cout << "The scheme (SLAU) is solved with a discrepancy ||R|| = " << I.getR(calcMatrVecMult(A, x), B) << std::endl;
			std::cout << "Estimation of the scheme error Z = " << I.getZ(u, x) << std::endl;
		}
	}
	*/

	void testModelTask() {
		std::cout << "Modul 7. Model task for estimating computational error\n";

		{
			this->prepareSerialDataForTest(5);
			ASSERT(node > 2);

			// Exact solution
			for (size_t i = 0; i < node; i++) {
				u[i] = 10 + 90 * x[i] * x[i];
			}

			Instrumental::printVec(x, "x");
			Instrumental::printVec(u, "u");
			print(h);

			double total = 12. / (h * h);

			a = total;
			c = 2 * total + 5.;
			b = total;

			print(a);
			print(b);
			print(c);

			vec Phi(node - 2);
			for (size_t i = 0; i < node - 2; i++) {
				Phi[i] = 450 * x[i + 1] * x[i + 1] - 2110.;
			}

			Instrumental::printVec(Phi, "Phi 1");

			size_t n = N - 2;
			vec A(n, a);
			vec B(n, b);
			vec C(n, c);
			pairs mu = std::make_pair(10., 100.);
			pairs kappa = std::make_pair(0., 0.);
			pairs gamma = std::make_pair(1., 1.);

			// Private solution
			SerialSweepMethod sweepMethod(A, C, B, Phi, kappa, mu, gamma);
			v = sweepMethod.run();

			I.setUV(u, v);

			ASSERT_EQUAL(v, u);

			// Numeric methods
			matr res = S.createMatr();

			Instrumental::printMatr(res, "res");

			Phi.insert(Phi.begin(), mu.first);
			Phi.insert(--Phi.end(), mu.second);

			Instrumental::printVec(Phi, "Phi 2");

			printf("The scheme (SLAU) is solved with a discrepancy ||R|| = %f\n", I.calcR(I.calcMatrVecMult(res, v), Phi));
			printf("Estimation of the scheme error Z = %f\n", I.calcZ());
		}
	}

	void execute() {
		TestRunner testRunner;
		str line = "-------------------------------";

		std::vector<std::function<void()>> tests = {
			[this]() {this->testModelTask();  }
			// [this]() { this->testEnteredData(); },
			// [this]() { this->test1(); },
			// [this]() { this->test2(); },
		};

		for (auto test : tests) {
			std::cout << line << std::endl;

			RUN_TEST(testRunner, test);

			std::cout << line << std::endl;
		}
	}
};

