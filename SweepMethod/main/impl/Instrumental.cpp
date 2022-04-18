#include "main/interfaces/Instrumental.h"
#include <test/common/TestRunner.h>

#include <ctime>


Instrumental::Instrumental(size_t n, TASK task_) {
    N = n;
    h = 1 / static_cast<double>(n);
    x = std::move(this->getGridNodes());

    u.assign(n, 0.);
    v.assign(n, 0.);

    task = task_;
}

vec Instrumental::getGridNodes() const {
    vec res(N);
    for (int i = 0; i < N; i++) {
        res[i] = (double)i * h;
    }

    return res;
}

std::tuple<size_t, double, vec, vec, vec, Task::TASK> Instrumental::getAllFields() const {
	return std::make_tuple(N, h,
                           x, u, v,
                           task
    );
}

void Instrumental::setAllFields(size_t n, double h_,
                                const vec& x_, const vec& v_, const vec& u_,
                                const Task::TASK& task_) {
    this->N = n;
    this->h = h_;
    this->x = x_;
    this->v = v_;
    this->u = u_;
    this->task = task_;
}


void Instrumental::printVec(const vec& a, const str& name) {
	if (a.size() < 30) {
		std::cout << name << ":\n";
		bool flag = true;
		for (double i : a) {
			if (!flag) {
				std::cout << ", ";
			} flag = false;
			std::cout << i;
		} std::cout << std::endl;
	}
}

void Instrumental::printMatr(const matr& a, const str& name) {
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

double Instrumental::calcR(const vec& x, const vec& b) const {
	double R = 0.;

	for (int i = 0; i < N; i++) {
		R = fmax(R, fabs(b[i] - x[i]));
	}

	return R;
}

vec Instrumental::calcMatrVecMult(const matr& A, const vec& b) {
	size_t n = A.size();
    size_t i, j;

	vec res(n);

	#pragma omp parallel for private(i, j) shared(res, n, A, b) default(none)
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            res[i] += A[i][j] * b[j];
        }
    }

	return res;
}

double Instrumental::calcZ() const {
	double Z = 0.;

	for (int i = 0; i < N; i++) {
		Z = fmax(Z, fabs(u[i] - v[i]));
	}

	return Z;
}

matr Instrumental::createThirdDiagMatr(double diag, double upDiag, double downDiag,
                                               double first, double second,
                                               double preLast, double last) {
    matr res(N, vec(N));

    #pragma omp parallel for shared(N, res, diag, upDiag, downDiag) default(none) if (N > 500)
    for (int i = 1; i < N; i++) {
        for (int j = 0; j < N; j++) {
            res[i][i] = diag;
            res[i][i - 1] = upDiag;
            res[i - 1][i] = downDiag;
        }
    }

    res[0][0] = first;  res[N - 1][N - 1] = last;
    res[0][1] = second; res[N - 1][N - 2] = preLast;

    return res;
}

bool Instrumental::compareDouble(double a, double b) {
    return std::fabs(a - b) <= EPS;
}

bool Instrumental::compareMatr(const matr& a, const matr& b) {
    size_t n = a.size();
    size_t i, j;

    #pragma omp parallel for private(i, j) shared(a, b, n) default(none)
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            Instrumental::compareDouble(a[i][j], b[i][j]);
        }
    }

    return true;
}

bool Instrumental::compareVec(const vec& a, const vec& b) {
    size_t n = a.size();
    size_t i;

    #pragma omp parallel for private(i) shared(a, b, n) default(none)
    for (i = 0; i < n; i++) {
        Instrumental::compareDouble(a[i], b[i]);
    }

    return true;
}

const vec &Instrumental::getU() const {
    return u;
}
