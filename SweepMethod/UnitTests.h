#pragma once

#include "TestRunner.h"
#include "Profiler.h"

void testPrepareData() {
	std::cout << "Проверка корректности введенных данных\n";
	int N, THREADNUM;

	{
		N = 8;
		THREADNUM = 4;
		ASSERT(!checkData(N, THREADNUM));
	}
	{
		N = 80;
		THREADNUM = 10;
		ASSERT(checkData(N, THREADNUM));
	}
	{
		N = 100;
		THREADNUM = 10;
		ASSERT(checkData(N, THREADNUM));
	}
}

void testSerialSweepMethod() {
	std::cout << "Проверка корректности работы последовательного метода прогонки\n";

	int N, THREADNUM, blockSize, classicSize, n;

	{
		std::tie(N, THREADNUM, blockSize, classicSize, n) = std::make_tuple(10, 4, 4, 6, 8);

		double h = 1 / static_cast<double>(N);
		vec x = getX(N);
		double a = 12 / pow(h, 2);
		double c = a * 2 + 5;

		vec F(n);
		vec A(n, a);
		vec B(n, a);
		vec C(n, c);
		pairs mu = std::make_pair(10, 100);
		pairs kappa = std::make_pair(0, 0);


		// частное решение
		vec v = serialSweepMethod(A, C, B, F, kappa, mu);

		// точное решение
		vec u(N);
		for (int i = 0; i < N; i++) {
			u[i] = 10 + 90 * pow(x[i], 2);
		}

		printVec(v, "v");
		printVec(u, "u");

		ASSERT_EQUAL(v, u);

		// численные методы
		matr res = createMatr(N, a, a, c);
		printMatr(res, "res");

		vec::iterator it = F.begin();
		F.insert(it, mu.first);
		it = F.end();
		F.insert(--it, mu.second);

		printVec(F, "fFull");
		std::cout << "Схема (СЛАУ) решена с невязкой ||R|| = " << getR(calcMatrVecMult(res, v), F) << std::endl;
		std::cout << "Оценка погрешности схемы Z = " << getZ(u, v) << std::endl;
	}
}

void test1() {
	int N, THREADNUM, blockSize, classicSize;
	std::tie(N, THREADNUM, blockSize, classicSize) = prepareData();
	int M1 = classicSize;
	matr A = createThirdDiagMatrI(N);
	vec  B = createVecN(N);

	std::cout << "----------------------" << std::endl;

	// Последовательная реализация прогонки
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

		vec x = serialSweepMethod(a, c, b, f, kappa, mu, gamma);

		printVec(x, "\nSerial x");

		std::cout << "Схема (СЛАУ) решена с невязкой ||R|| = " << getR(calcMatrVecMult(A, x), B) << std::endl;
	}

	std::cout << "----------------------" << std::endl;

	// Параллельная реализация прогонки
	{
		LOG_DURATION("ParallelSweepMethod");

		vec x = ParallelSweepMethod(N, THREADNUM, blockSize, M1);

		printVec(x, "\nParallel x");

		std::cout << "Схема (СЛАУ) решена с невязкой ||R|| = " << getR(calcMatrVecMult(A, x), B) << std::endl;
	}

	std::cout << "----------------------" << std::endl;
}

void test2() {
	int N, THREADNUM, blockSize, classicSize;
	std::tie(N, THREADNUM, blockSize, classicSize) = prepareData();
	int M1 = classicSize;

	matr A = createThirdDiagMatrI(N);
	vec  B = createVecN(N);

	std::cout << "----------------------" << std::endl;

	// Последовательная реализация прогонки (с точным решением)
	{
		LOG_DURATION("SerialSweepMethod");

		vec xi(N);
		double h = 1 / N;
		for (int i = 0; i < N; i++) {
			xi[i] = i * h;
		}

		/*
		std::cout << "Схема (СЛАУ) решена с невязкой ||R|| = " << getR(calcMatrVecMult(A, x), B) << std::endl;
		std::cout << "Оценка погрешности схемы Z = " << getZ(u, x) << std::endl;
		*/
	}
}

void callAllTests() {
	TestRunner testRunner;

	/*
	RUN_TEST(testRunner, test1);
	RUN_TEST(testRunner, test2);
	*/

	// RUN_TEST(testRunner, testPrepareData);
	RUN_TEST(testRunner, testSerialSweepMethod);
}