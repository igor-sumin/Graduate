#pragma once

#include <functional>

#include "TestRunner.h"
#include "Profiler.h"

#include "main/interfaces/tools/Instrumental.h"
#include "main/interfaces/serial/SerialSweepMethod.h"
#include "main/interfaces/parallel/ParallelSweepMethod.h"

class UnitTests {
private:
	int N, threadNum, blockSize, classicSize;
	matr A;
	vec B;
	Instrumental I;

	void prepareDataForTest() {
		A = I.createThirdDiagMatrI();
		B = I.createVecN();

		tie(N, threadNum, blockSize, classicSize) = I.getFields();
	}

public:
	UnitTests() : N(10), threadNum(4), blockSize(4), classicSize(6) {}

	void testEnteredData() {
		std::cout << "Checking the correctness of the entered data\n";

		{
			Instrumental i(8, 4, -1, -1);
			ASSERT(!i.checkData());
		}
		{
			Instrumental i(80, 10, -1, -1);
			ASSERT(!i.checkData());
		}
		{
			Instrumental i(100, 10, -1, -1);
			ASSERT(!i.checkData());
		}
	}

	void testSerialSweepMethod() {
		std::cout << "Checking the correctness of the sequential sweep method\n";

		int N, threadNum, blockSize, classicSize, n;
		{
			Instrumental instrumentalTest(10, 4, 4, 6);
			std::tie(N, threadNum, blockSize, classicSize) = instrumentalTest.getFields();

			double h = 1 / static_cast<double>(N);
			vec x = instrumentalTest.getX();
			double a = 12 / pow(h, 2);
			double c = a * 2 + 5;

			n = N - 2;
			vec F(n);
			vec A(n, a);
			vec B(n, a);
			vec C(n, c);
			pairs mu = std::make_pair(10, 100);
			pairs kappa = std::make_pair(0, 0);
			pairs gamma = std::make_pair(1, 1);

			// Private solution
			SerialSweepMethod ssm(A, C, B, F, kappa, mu, gamma);
			vec v = ssm.run();

			// Exact solution
			vec u(N);
			for (int i = 0; i < N; i++) {
				u[i] = 10 + 90 * pow(x[i], 2);
			}

			Instrumental::printVec(v, "v");
			Instrumental::printVec(u, "u");

			ASSERT_EQUAL(v, u);

			// Numeric methods
			matr res = instrumentalTest.createMatr(a, a, c);
			Instrumental::printMatr(res, "res");

			vec::iterator it = F.begin();
			F.insert(it, mu.first);
			it = F.end();
			F.insert(--it, mu.second);

			Instrumental::printVec(F, "fFull");

			std::cout << "The scheme (SLAU) is solved with a discrepancy ||R|| = "
				<< instrumentalTest.getR(instrumentalTest.calcMatrVecMult(res, v), F) << std::endl;

			std::cout << "Estimation of the scheme error Z = "
				<< instrumentalTest.getZ(u, v) << std::endl;
		}
	}

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
			*/
		}
	}

	void execute() {
		TestRunner testRunner;
		str line = "-------------------------------";

		std::vector<std::function<void()>> tests = {
			[this]() { this->testEnteredData(); },
			[this]() { this->test1(); },
			[this]() { this->test2(); },
			[this]() { this->testSerialSweepMethod(); }
		};

		for (auto test : tests) {
			std::cout << line << std::endl;

			RUN_TEST(testRunner, test);

			std::cout << line << std::endl;
		}
	}
};