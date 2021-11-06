#pragma once

void transformMatrA(matr& A, int N, int THREADNUM, int blockSize) {
	for (int iter = 0; iter < THREADNUM; iter++) {
		// сверху-вниз
		for (int i = iter * blockSize; i < (iter + 1) * blockSize - 1; i++) {
			double k = A[i + 1][i] / A[i][i];
			for (int j = 0; j < i + 1; j++) {
				A[i + 1][j] -= k * A[i][j];
			}
		}

		// снизу-вверх
		for (int i = (iter + 1) * blockSize - 2; i >= iter * blockSize; i--) {
			double k = A[i][i + 1] / A[i + 1][i + 1];
			for (int j = i; j < N - 1; j++) {
				A[i][j + 1] -= k * A[i + 1][j + 1];
			}
		}
	}
}

std::pair<matr, vec> collectInterferElem(matr& A, vec& b, int N, int THREADNUM, int blockSize, int M1) {
	matr R(M1, vec(M1, 0.));
	vec  Y1(M1);

	// Находим вектор неизвестных
	for (int i = blockSize, k = 1; i < N; i += blockSize, k += 2) {
		Y1[k - 1] = b[i - 1];
		Y1[k] = b[i];
	}

	// 1. предобработка (крайние части)
	// -> верхняя
	int Nb = blockSize - 1;
	R[0][0] = A[Nb][Nb];
	R[0][1] = A[Nb][blockSize];
	R[1][0] = A[blockSize][Nb];
	R[2][1] = A[2 * blockSize - 1][Nb];

	// -> нижняя
	Nb = N - blockSize;
	R[M1 - 1][M1 - 1] = A[Nb][Nb];
	R[M1 - 1][M1 - 2] = A[Nb][Nb - 1];
	R[M1 - 2][M1 - 1] = A[Nb - 1][Nb];
	R[M1 - 3][M1 - 2] = A[N - blockSize * 2][Nb];

	// 2. постобработка (внутренняя часть)
	for (int i = blockSize, k = 1; i < N - blockSize; i += blockSize, k += 2) {
		for (int j = blockSize, l = 1; j < N - blockSize; j += blockSize, l += 2) {
			double a = A[i][j];
			double b = A[i + blockSize - 1][j];
			double c = A[i][j + blockSize - 1];
			double d = A[i + blockSize - 1][j + blockSize - 1];

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
	for (int i = 0; i < M1; i += 2) {
		std::swap(R[i][i], R[i][i + 1]);
		std::swap(R[i + 1][i], R[i + 1][i + 1]);
		std::swap(Y1[i], Y1[i + 1]);
	}


	return std::make_pair(R, Y1);
}

vec ClassicSweepMethod(const matr& R, const vec& Y1, int M1) {
	vec a(M1 - 2),
		c(M1 - 2),
		b(M1 - 2),
		f(M1 - 2);

	pairs mu    = std::make_pair(Y1[0], Y1[M1 - 1]);
	pairs kappa = std::make_pair(-R[0][1], -R[M1 - 1][M1 - 2]);
	pairs gamma = std::make_pair(R[0][0] , R[M1 - 1][M1 - 1]);


	#pragma omp parallel sections
	{
		#pragma omp section
		{
			for (int i = 1; i < M1 - 1; i++) {
				f[i - 1] = -Y1[i];
				a[i - 1] = R[i][i - 1];
			}
		}

		#pragma omp section
		{
			for (int i = 1; i < M1 - 1; i++) {
				b[i - 1] = R[i][i + 1];
				c[i - 1] = -R[i][i];
			}
		}
	}

	return serialSweepMethod(a, c, b, f, kappa, mu, gamma);
}

// Получаем Y2 из векторов b, Y1
vec getY2(vec& b, vec& Y1) {
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

vec getX2(const matr& A, vec& Y2, int N, int blockSize, int M1) {
	vec X2(N - M1);

	// 1. предобработка (крайние части)
	// -> верхняя + нижняя
	int last = N - blockSize - 1;
	for (int i = 0; i < blockSize - 1; i++) {
		int j = N - i - 1;

		// находим коэф-ты
		Y2[i] -= A[i][blockSize];
		Y2[j - M1] -= A[j][last];

		// находим вектор неизвестных
		X2[i] = Y2[i] / A[i][i];
		X2[j - M1] = Y2[j - M1] / A[j][j];
	}

	// 2. постобработка (внутреняя часть)
	for (int k = blockSize + 1, l = 1; k < N - blockSize - 1; k += blockSize, l++) {
		for (int i = 0; i < blockSize - 2; i++) {
			// находим коэф-ты
			Y2[i + k - 2 * l] -= (A[i + k][k - 2] + A[i + k][k + blockSize - 1]);

			// находим вектор неизвестных
			X2[i + k - 2 * l] = Y2[i + k - 2 * l] / A[i + k][i + k];
		}
	}

	return X2;
}

vec result(const vec& X1, const vec& X2, int N, int blockSize, int M1) {
	vec X(N);

	// 1. предобработка (крайние части)
	// -> верхняя + нижняя
	for (int i = 0; i < blockSize - 1; i++) {
		int j = N - i - 1;

		X[i] = X2[i];
		X[j] = X2[j - M1];
	}

	int i = 0,
		j = blockSize - 1,
		k = blockSize - 1;

	// 2. постобработка (внутреняя часть)
	while (k < N - blockSize) {
		for (int jiter = 0; jiter < 2; jiter++) {
			X[k++] = X1[i++];
		}

		for (int kiter = 0; kiter < blockSize - 2; kiter++) {
			X[k++] = X2[j++];
		}
	}

	return X;
}

// Вызов всех методов
vec ParallelSweepMethod(int N, int THREADNUM, int blockSize, int M1) {
	// 0. Подготовка данных
	matr A = createThirdDiagMatrI(N);
	vec  b = createVecN(N);

	printMatr(A, "A1");
	printVec(b, "b");

	// 1. Приведение матрицы к столбцовым видам
	transformMatrA(A, N, THREADNUM, blockSize);

	printMatr(A, "A2");

	// 2. Формируем матрицу R для последовательной прогонки вместе с вектором неизвестных Y1
	matr R; vec Y1;
	std::tie(R, Y1) = collectInterferElem(A, b, N, THREADNUM, blockSize, M1);

	printMatr(R, "R2");
	printVec(Y1, "Y1");

	// 3. Решаем СЛАУ (R*X1 = Y1) классической прогонкой
	vec X1 = ClassicSweepMethod(R, Y1, M1);
	printVec(X1, "X1");

	// 4. Получаем вектор неизвестных Y2
	vec Y2 = getY2(b, Y1);
	printVec(Y2, "Y2");

	// 5. Находим вектор неизвестных X2
	vec X2 = getX2(A, Y2, N, blockSize, M1);
	printVec(X2, "X2");

	// 6. Группируем вектора неизвестных X1, X2 в X
	return result(X1, X2, N, blockSize, M1);
}