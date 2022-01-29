#include <test/common/Profiler.h>
#include "main/interfaces/parallel/ParallelSweepMethod.h"


std::tuple<size_t, size_t, size_t, size_t, matr, vec, vec> ParallelSweepMethod::getAllFields() const {
    return std::make_tuple(N, threadNum, blockSize, interSize, A, b, y);
}

void ParallelSweepMethod::setAllFields(size_t N, size_t threadNum, size_t blockSize, size_t interSize, matr A, vec b, vec y) {
    this->N = N;
    this->threadNum = threadNum;
    this->blockSize = blockSize;
    this->interSize = interSize;
    this->A = std::move(A);
    this->b = std::move(b);
    this->y = std::move(y);
}

void ParallelSweepMethod::transformation() {
    size_t i, j, k;
    double coef;
    int iter;

    #pragma omp parallel private(i, j, k, coef, iter) shared(A, b, blockSize, N) default(none)
    {
        iter = omp_get_thread_num();

        // top-down
        for (i = iter * blockSize + 1; i < (iter + 1) * blockSize; i++) {
            coef = A[i][i - 1] / A[i - 1][i - 1];
            for (j = 0; j < i + 1; j++) {
                A[i][j] -= coef * A[i - 1][j];
            }
            b[i] -= coef * b[i - 1];
        }

        // bottom-up
        for (k = (iter + 1) * blockSize - 1; k > iter * blockSize; k--) {
            i = k - 1;
            coef = A[i][i + 1] / A[i + 1][i + 1];
            for (j = i; j < N - 1; j++) {
                A[i][j + 1] -= coef * A[i + 1][j + 1];
            }
            b[i] -= coef * b[i + 1];
        }
    }
}

void ParallelSweepMethod::preULR(matr& R) {
    size_t bS1 = blockSize - 1;
    R[0][0] = A[bS1][bS1];
    R[0][1] = A[bS1][blockSize];
    R[1][0] = A[blockSize][bS1];
    R[2][1] = A[2 * blockSize - 1][bS1];
}

void ParallelSweepMethod::preLRR(matr& R) {
    size_t bSN = N - blockSize;
    R[interSize - 1][interSize - 1] = A[bSN][bSN];
    R[interSize - 1][interSize - 2] = A[bSN][bSN - 1];
    R[interSize - 2][interSize - 1] = A[bSN - 1][bSN];
    R[interSize - 3][interSize - 2] = A[N - blockSize * 2][bSN];
}

std::pair<matr, vec> ParallelSweepMethod::collectInterferElem() {
    matr R(interSize, vec(interSize, 0.));
    vec partB(interSize);

    size_t k = 1, l = 1;
    size_t iter, jter;
    size_t i, j, s;
    double a1, a2, a3, a4;

    // 1. pre-processing (extreme part)
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

    // 2. post-processing (internal part)
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

	// ordering the coefficient
    #pragma omp parallel for private(i) shared(R, partB, interSize) default(none)
    for (i = 0; i < interSize; i += 2) {
        std::swap(R[i][i], R[i][i + 1]);
        std::swap(R[i + 1][i], R[i + 1][i + 1]);
    }

	return std::make_pair(R, partB);
}

vec ParallelSweepMethod::collectPartY(const matr& R, const vec& partB) {
    vec a_(interSize - 2),
        c_(interSize - 2),
        b_(interSize - 2),
        phi(interSize - 2);

    size_t i;
    pairs kappa = std::make_pair(-R[0][1], -R[interSize - 1][interSize - 2]);
    pairs mu = std::make_pair(partB[0], partB[interSize - 1]);
    pairs gamma = std::make_pair(R[0][0], R[interSize - 1][interSize - 1]);

    #pragma omp parallel for private(i) shared(a_, c_, b_, phi, partB, R, interSize) default(none)
    for (i = 1; i < interSize - 1; i++) {
        a_[i - 1] = R[i][i - 1];
        b_[i - 1] = R[i][i + 1];
        c_[i - 1] = -R[i][i];
        phi[i - 1] = -partB[i];
    }

    SerialSweepMethod ssm(a_, c_, b_, phi, kappa, mu, gamma);
    vec partY = ssm.run();

    #pragma omp parallel for private(i) shared(partY, interSize) default(none)
    for (i = 0; i < interSize - 1; i += 2) {
        std::swap(partY[i], partY[i + 1]);
    }

    return partY;
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
	vec X2(N - interSize);

    // 1. preprocessing (extreme parts)
    // -> upper + lower
	int last = N - blockSize - 1;
	for (int i = 0; i < blockSize - 1; i++) {
		int j = N - i - 1;

		// finding coefficients
		Y2[i] -= A[i][blockSize];
		Y2[j - interSize] -= A[j][last];

		// finding vector of unknowns
		X2[i] = Y2[i] / A[i][i];
		X2[j - interSize] = Y2[j - interSize] / A[j][j];
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
		X[j] = X2[j - interSize];
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
    transformation();

	printMatr(A, "A2");

    // 2. From the matrix R for a serial sweep method along with the vector of unknowns partB
	matr R; vec Y1;
	std::tie(R, Y1) = collectInterferElem();

	printMatr(R, "R2");
	printVec(Y1, "partB");

	// 3. Solve SLAU (R * X1 = partB) by the classical sweep method
	vec X1 = collectPartY(R, Y1);
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