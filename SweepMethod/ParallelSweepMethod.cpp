#include "ParallelSweepMethod.h"


void ParallelSweepMethod::transformMatrA() {
	for (int iter = 0; iter < THREAD_NUM; iter++) {
		// сверху-вниз
		for (int i = iter * BLOCK_SIZE; i < (iter + 1) * BLOCK_SIZE - 1; i++) {
			double k = A[i + 1][i] / A[i][i];
			for (int j = 0; j < i + 1; j++) {
				A[i + 1][j] -= k * A[i][j];
			}
		}

		// снизу-вверх
		for (int i = (iter + 1) * BLOCK_SIZE - 2; i >= iter * BLOCK_SIZE; i--) {
			double k = A[i][i + 1] / A[i + 1][i + 1];
			for (int j = i; j < N - 1; j++) {
				A[i][j + 1] -= k * A[i + 1][j + 1];
			}
		}
	}
}

std::pair<matr, vec> ParallelSweepMethod::collectInterferElem() {
	matr R(CLASSIC_SIZE, vec(CLASSIC_SIZE, 0.));
	vec  Y1(CLASSIC_SIZE);

	// Находим вектор неизвестных
	for (int i = BLOCK_SIZE, k = 1; i < N; i += BLOCK_SIZE, k += 2) {
		Y1[k - 1] = b[i - 1];
		Y1[k] = b[i];
	}

	// 1. предобработка (крайние части)
	// -> верхняя
	int Nb = BLOCK_SIZE - 1;
	R[0][0] = A[Nb][Nb];
	R[0][1] = A[Nb][BLOCK_SIZE];
	R[1][0] = A[BLOCK_SIZE][Nb];
	R[2][1] = A[2 * BLOCK_SIZE - 1][Nb];

	// -> нижняя
	Nb = N - BLOCK_SIZE;
	R[CLASSIC_SIZE - 1][CLASSIC_SIZE - 1] = A[Nb][Nb];
	R[CLASSIC_SIZE - 1][CLASSIC_SIZE - 2] = A[Nb][Nb - 1];
	R[CLASSIC_SIZE - 2][CLASSIC_SIZE - 1] = A[Nb - 1][Nb];
	R[CLASSIC_SIZE - 3][CLASSIC_SIZE - 2] = A[N - BLOCK_SIZE * 2][Nb];

	// 2. постобработка (внутренняя часть)
	for (int i = BLOCK_SIZE, k = 1; i < N - BLOCK_SIZE; i += BLOCK_SIZE, k += 2) {
		for (int j = BLOCK_SIZE, l = 1; j < N - BLOCK_SIZE; j += BLOCK_SIZE, l += 2) {
			double a = A[i][j];
			double b = A[i + BLOCK_SIZE - 1][j];
			double c = A[i][j + BLOCK_SIZE - 1];
			double d = A[i + BLOCK_SIZE - 1][j + BLOCK_SIZE - 1];

			if (a != 0 && d != 0) {
				R[k][l] = a;
				R[k + 1][l + 1] = d;
			}
			else if (a != 0) {
				R[k][l - 1] = a;
				R[k + 1][l] = b;
			}
			else if (d != 0) {
				R[k][l + 1] = c;
				R[k + 1][l + 2] = d;
			}
		}
	}

	printMatr(R, "R1");
	printVec(Y1, "Y1");

	// упорядочиваем коэф-ты
	for (int i = 0; i < CLASSIC_SIZE; i += 2) {
		std::swap(R[i][i], R[i][i + 1]);
		std::swap(R[i + 1][i], R[i + 1][i + 1]);
		std::swap(Y1[i], Y1[i + 1]);
	}


	return std::make_pair(R, Y1);
}

vec ParallelSweepMethod::ClassicSweepMethod(const matr& R, const vec& Y1) {
	vec a(CLASSIC_SIZE - 2),
		c(CLASSIC_SIZE - 2),
		b(CLASSIC_SIZE - 2),
		f(CLASSIC_SIZE - 2);

	pairs mu = std::make_pair(Y1[0], Y1[CLASSIC_SIZE - 1]);
	pairs kappa = std::make_pair(-R[0][1], -R[CLASSIC_SIZE - 1][CLASSIC_SIZE - 2]);
	pairs gamma = std::make_pair(R[0][0], R[CLASSIC_SIZE - 1][CLASSIC_SIZE - 1]);


	#pragma omp parallel sections
	{
		#pragma omp section
		{
			for (int i = 1; i < CLASSIC_SIZE - 1; i++) {
				f[i - 1] = -Y1[i];
				a[i - 1] = R[i][i - 1];
			}
		}

		#pragma omp section
		{
			for (int i = 1; i < CLASSIC_SIZE - 1; i++) {
				b[i - 1] = R[i][i + 1];
				c[i - 1] = -R[i][i];
			}
		}
	}

	SerialSweepMethod ssm(a, c, b, f, kappa, mu, gamma);
	return ssm.run();
}

// Получаем Y2 из векторов b, Y1
vec ParallelSweepMethod::getY2(vec& Y1) {
	vec Y2;

	std::copy_if(b.begin(), b.end(),
		std::back_inserter(Y2),
		std::bind(
			std::equal_to<>(),
			std::bind(
				std::find<std::vector<double>::iterator, double>,
				Y1.begin(), Y1.end(), std::placeholders::_1
			), Y1.end()
		)
	);

	return Y2;
}

vec ParallelSweepMethod::getX2(vec& Y2) {
	vec X2(N - CLASSIC_SIZE);

	// 1. предобработка (крайние части)
	// -> верхняя + нижняя
	int last = N - BLOCK_SIZE - 1;
	for (int i = 0; i < BLOCK_SIZE - 1; i++) {
		int j = N - i - 1;

		// находим коэф-ты
		Y2[i] -= A[i][BLOCK_SIZE];
		Y2[j - CLASSIC_SIZE] -= A[j][last];

		// находим вектор неизвестных
		X2[i] = Y2[i] / A[i][i];
		X2[j - CLASSIC_SIZE] = Y2[j - CLASSIC_SIZE] / A[j][j];
	}

	// 2. постобработка (внутреняя часть)
	for (int k = BLOCK_SIZE + 1, l = 1; k < N - BLOCK_SIZE - 1; k += BLOCK_SIZE, l++) {
		for (int i = 0; i < BLOCK_SIZE - 2; i++) {
			// находим коэф-ты
			Y2[i + k - 2 * l] -= (A[i + k][k - 2] + A[i + k][k + BLOCK_SIZE - 1]);

			// находим вектор неизвестных
			X2[i + k - 2 * l] = Y2[i + k - 2 * l] / A[i + k][i + k];
		}
	}

	return X2;
}

vec ParallelSweepMethod::result(const vec& X1, const vec& X2) {
	vec X(N);

	// 1. предобработка (крайние части)
	// -> верхняя + нижняя
	for (int i = 0; i < BLOCK_SIZE - 1; i++) {
		int j = N - i - 1;

		X[i] = X2[i];
		X[j] = X2[j - CLASSIC_SIZE];
	}

	int i = 0,
		j = BLOCK_SIZE - 1,
		k = BLOCK_SIZE - 1;

	// 2. постобработка (внутреняя часть)
	while (k < N - BLOCK_SIZE) {
		for (int jiter = 0; jiter < 2; jiter++) {
			X[k++] = X1[i++];
		}

		for (int kiter = 0; kiter < BLOCK_SIZE - 2; kiter++) {
			X[k++] = X2[j++];
		}
	}

	return X;
}

// Вызов всех методов
vec ParallelSweepMethod::run() {
	printMatr(A, "A1");
	printVec(b, "b");

	// 1. Приведение матрицы к столбцовым видам
	transformMatrA();

	printMatr(A, "A2");

	// 2. Формируем матрицу R для последовательной прогонки вместе с вектором неизвестных Y1
	matr R; vec Y1;
	std::tie(R, Y1) = collectInterferElem();

	printMatr(R, "R2");
	printVec(Y1, "Y1");

	// 3. Решаем СЛАУ (R*X1 = Y1) классической прогонкой
	vec X1 = ClassicSweepMethod(R, Y1);
	printVec(X1, "X1");

	// 4. Получаем вектор неизвестных Y2
	vec Y2 = getY2(Y1);
	printVec(Y2, "Y2");

	// 5. Находим вектор неизвестных X2
	vec X2 = getX2(Y2);
	printVec(X2, "X2");

	// 6. Группируем вектора неизвестных X1, X2 в X
	return result(X1, X2);
}