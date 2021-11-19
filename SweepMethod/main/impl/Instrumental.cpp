#include "main/interfaces/Instrumental.h"

#include <ctime>


void Instrumental::prepareData() {
}

bool Instrumental::checkData() const {
	return true;
}

void Instrumental::setN(size_t n) {
	N = n;
}

void Instrumental::setUV(vec& u_, vec& v_) {
	Instrumental::printVec(u_, "u_");
	Instrumental::printVec(v_, "v_");

	u = u_;
	v = v_;
}

std::tuple<size_t, size_t, vec, vec> Instrumental::getAllFields() const {
	return std::make_tuple(N, node, u, v);
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

	for (int i = 0; i < node; i++) {
		R = fmax(R, fabs(b[i] - x[i]));
	}

	return R;
}

vec Instrumental::calcMatrVecMult(const matr& A, const vec& b) {
	int n = A.size();
	vec res(n);

	#pragma omp parallel shared(res, n, A, b) default(none)
	{
		#pragma omp for
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				res[i] += A[i][j] * b[j];
	}

	return res;
}

double Instrumental::calcZ() const {
	double Z = 0.;

	for (int i = 0; i < node; i++) {
		Z = fmax(Z, fabs(u[i] - v[i]));
	}

	return Z;
}