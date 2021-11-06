#include <tuple>
#include <iostream>

#include "Instrumental.h"

std::tuple<int, int, int, int> Instrumental::prepareData() {
	do {
		std::cout << "Введите размерность (N) = ";
		std::cin >> N;

		std::cout << "Введите количество вычислительных узлов (THREADNUM) = ";
		std::cin >> THREAD_NUM;

	} while (!checkData());

	BLOCK_SIZE   = N / THREAD_NUM;
	CLASSIC_SIZE = (THREAD_NUM - 1) * 2;

	return std::make_tuple(N, THREAD_NUM, BLOCK_SIZE, CLASSIC_SIZE);
}

bool Instrumental::checkData() {
	if (N < 7) {
		std::cout << "Размерность " << N << " слишком мала\n";
		std::cout << "Параллельные вычисления для нее не эффективны\n";
		std::cout << "Введите бОльшую размерность N\n\n";

		return false;
	}

	if ((N / THREAD_NUM) < 3) {
		std::cout << "Алгоритм нерабочий для пропорций размерности с количеством вычислительных узлов\n";
		std::cout << "Введите бОльшую размерность, либо уменьшите количество вычислительных узлов\n\n";

		return false;
	}

	vec div = findDivisors(N);
	if (isPrime(N) || std::find(div.begin(), div.end(), THREAD_NUM) == div.end()) {
		std::cout << "Невозможно разбить на одинаковые блоки размерность " << N << "\n";
		std::cout << "Введите корректные значения N, THREADNUM\n\n";

		return false;
	}

	return true;
}

bool Instrumental::isPrime(int num) const {
	bool ok = true;
	for (int i = 2; i <= num / 2; ++i) {
		if (num % i == 0) {
			ok = false;
			break;
		}
	}

	return ok;
}

vec Instrumental::findDivisors(int num) {
	vec res;

	for (int i = 2; i <= sqrt(num); i++) {
		if (num % i == 0) {
			if (num / i != i) {
				res.push_back(num / i);
			}

			res.push_back(i);
		}
	}

	return res;
}

void Instrumental::printVec(const vec& a, const std::string& name) {
	if (a.size() < 30) {
		std::cout << name << ":\n";
		bool flag = true;
		for (int i = 0; i < a.size(); i++) {
			if (!flag) {
				std::cout << ", ";
			} flag = false;
			std::cout << a[i];
		} std::cout << std::endl;
	}
}

void Instrumental::printMatr(const matr& a, const std::string& name) {
	if (a.size() < 30) {
		bool first = true;

		std::cout << name << ":\n";
		for (int i = 0; i < a.size(); i++) {
			first = true;
			for (int j = 0; j < a[0].size(); j++) {
				if (!first) {
					std::cout << ", ";
				} first = false;

				printf("%8.3f", a[i][j]);
			} printf("\n");
		} std::cout << "\n";
	}
}

vec Instrumental::createVecN() {
	vec a(N);
	std::iota(a.begin(), a.end(), 0);

	return a;
}

vec Instrumental::createVecRand() {
	std::mt19937 gen;
	gen.seed(static_cast<unsigned int>(time(0)));
	vec a(N);

	#pragma omp parallel for
	for (int i = 0; i < N; i++)
		a[i] = gen() % 100;

	return a;
}

matr Instrumental::createThirdDiagMatrI() {
	matr a(N, vec(N));

	#pragma omp parallel for if (N > 500)
	for (int i = 1; i < N; i++) {
		for (int j = 0; j < N; j++) {
			a[i][i] = 3.;
			a[i][i - 1] = 1.;
			a[i - 1][i] = 2.;
		}
	}

	a[0][0] = 1.; a[N - 1][N - 1] = 1.;
	a[0][1] = -0.5; a[N - 1][N - 2] = -0.5;

	return a;
}

matr Instrumental::createThirdDiagMatrRand() {
	std::mt19937 gen;
	gen.seed(static_cast<unsigned int>(time(0)));
	matr a(N, vec(N));

	a[0][0] = gen() % 100;
	#pragma omp parallel for if (n > 500)
	for (int i = 1; i < N; i++) {
		for (int j = 0; j < N; j++) {
			a[i][i] = gen() % 100;
			a[i][i - 1] = gen() % 100;
			a[i - 1][i] = gen() % 100;
		}
	}

	return a;
}

matr Instrumental::createMatr(double a, double b, double c) {
	matr res(N, vec(N));

	#pragma omp parallel for if (N > 500)
	for (int i = 1; i < N; i++) {
		for (int j = 0; j < N; j++) {
			res[i][i] = c;
			res[i][i - 1] = b;
			res[i - 1][i] = a;
		}
	}

	res[0][0] = 1.; res[N - 1][N - 1] = 1.;
	res[0][1] = 0.; res[N - 1][N - 2] = 0.;

	return res;
}

vec Instrumental::calcMatrVecMult(const matr& A, const vec& b) {
	int n = A.size();
	vec res(n);

#pragma omp parallel shared(res, n)
	{
#pragma omp for
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				res[i] += A[i][j] * b[j];
	}

	return res;
}