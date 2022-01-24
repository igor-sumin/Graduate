#include "main/interfaces/parallel/ParallelSweepMethod.h"


std::tuple<size_t, size_t, size_t, size_t, matr, vec> ParallelSweepMethod::getAllFields() const {
    return std::make_tuple(N, threadNum, blockSize, classicSize, A, b);
}

void ParallelSweepMethod::setAllFields(size_t N, size_t threadNum, size_t blockSize, size_t classicSize, matr A, vec b) {
    this->N = N;
    this->threadNum = threadNum;
    this->blockSize = blockSize;
    this->classicSize = classicSize;
    this->A = std::move(A);
    this->b = std::move(b);
}

void ParallelSweepMethod::transformMatrA() {
	for (int iter = 0; iter < threadNum; iter++) {
		// top-down
		for (int i = iter * blockSize; i < (iter + 1) * blockSize - 1; i++) {
			double k = A[i + 1][i] / A[i][i];
			for (int j = 0; j < i + 1; j++) {
				A[i + 1][j] -= k * A[i][j];
			}
		}

		// bottom-up
		for (int i = (iter + 1) * blockSize - 2; i >= iter * blockSize; i--) {
			double k = A[i][i + 1] / A[i + 1][i + 1];
			for (int j = i; j < N - 1; j++) {
				A[i][j + 1] -= k * A[i + 1][j + 1];
			}
		}
	}
}

std::pair<matr, vec> ParallelSweepMethod::collectInterferElem() {
	matr R(classicSize, vec(classicSize, 0.));
	vec  Y1(classicSize);

	// finding vector of unknowns
	for (int i = blockSize, k = 1; i < N; i += blockSize, k += 2) {
		Y1[k - 1] = b[i - 1];
		Y1[k] = b[i];
	}

    // 1. preprocessing (extreme parts)
    // -> upper
	int Nb = blockSize - 1;
	R[0][0] = A[Nb][Nb];
	R[0][1] = A[Nb][blockSize];
	R[1][0] = A[blockSize][Nb];
	R[2][1] = A[2 * blockSize - 1][Nb];

	// -> lower
	Nb = N - blockSize;
	R[classicSize - 1][classicSize - 1] = A[Nb][Nb];
	R[classicSize - 1][classicSize - 2] = A[Nb][Nb - 1];
	R[classicSize - 2][classicSize - 1] = A[Nb - 1][Nb];
	R[classicSize - 3][classicSize - 2] = A[N - blockSize * 2][Nb];

	// 2. post-processing (internal part)
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

	// ordering the coefficient
	for (int i = 0; i < classicSize; i += 2) {
		std::swap(R[i][i], R[i][i + 1]);
		std::swap(R[i + 1][i], R[i + 1][i + 1]);
		std::swap(Y1[i], Y1[i + 1]);
	}


	return std::make_pair(R, Y1);
}

vec ParallelSweepMethod::ClassicSweepMethod(const matr& R, const vec& Y1) {
	vec a(classicSize - 2),
		c(classicSize - 2),
		b(classicSize - 2),
		f(classicSize - 2);

	pairs mu = std::make_pair(Y1[0], Y1[classicSize - 1]);
	pairs kappa = std::make_pair(-R[0][1], -R[classicSize - 1][classicSize - 2]);
	pairs gamma = std::make_pair(R[0][0], R[classicSize - 1][classicSize - 1]);


	#pragma omp parallel sections
	{
		#pragma omp section
		{
			for (int i = 1; i < classicSize - 1; i++) {
				f[i - 1] = -Y1[i];
				a[i - 1] = R[i][i - 1];
			}
		}

		#pragma omp section
		{
			for (int i = 1; i < classicSize - 1; i++) {
				b[i - 1] = R[i][i + 1];
				c[i - 1] = -R[i][i];
			}
		}
	}

	SerialSweepMethod ssm(a, c, b, f, kappa, mu, gamma);
	return ssm.run();
}

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
	vec X2(N - classicSize);

    // 1. preprocessing (extreme parts)
    // -> upper + lower
	int last = N - blockSize - 1;
	for (int i = 0; i < blockSize - 1; i++) {
		int j = N - i - 1;

		// finding coefficients
		Y2[i] -= A[i][blockSize];
		Y2[j - classicSize] -= A[j][last];

		// finding vector of unknowns
		X2[i] = Y2[i] / A[i][i];
		X2[j - classicSize] = Y2[j - classicSize] / A[j][j];
	}

	// 2. post-processing (internal part)
	for (int k = blockSize + 1, l = 1; k < N - blockSize - 1; k += blockSize, l++) {
		for (int i = 0; i < blockSize - 2; i++) {
			// finding coefficients
			Y2[i + k - 2 * l] -= (A[i + k][k - 2] + A[i + k][k + blockSize - 1]);

			// finding vector of unknowns
			X2[i + k - 2 * l] = Y2[i + k - 2 * l] / A[i + k][i + k];
		}
	}

	return X2;
}

vec ParallelSweepMethod::result(const vec& X1, const vec& X2) {
	vec X(N);

    // 1. preprocessing (extreme parts)
    // -> upper + lower
	for (int i = 0; i < blockSize - 1; i++) {
		int j = N - i - 1;

		X[i] = X2[i];
		X[j] = X2[j - classicSize];
	}

	int i = 0,
		j = blockSize - 1,
		k = blockSize - 1;

    // 2. post-processing (internal part)
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

vec ParallelSweepMethod::run() {
	printMatr(A, "A1");
	printVec(b, "b");

    // 1. Reduction of the matrix to columnar views
	transformMatrA();

	printMatr(A, "A2");

    // 2. From the matrix R for a serial sweep method along with the vector of unknowns Y1
	matr R; vec Y1;
	std::tie(R, Y1) = collectInterferElem();

	printMatr(R, "R2");
	printVec(Y1, "Y1");

	// 3. Solve SLAU (R * X1 = Y1) by the classical sweep method
	vec X1 = ClassicSweepMethod(R, Y1);
	printVec(X1, "X1");

	// 4. Get the vector of unknowns Y2
	vec Y2 = getY2(Y1);
	printVec(Y2, "Y2");

	// 5. Find the vector of unknowns X2
	vec X2 = getX2(Y2);
	printVec(X2, "X2");

	// 6. Group vectors of unknowns X1, X2 into X
	return result(X1, X2);
}