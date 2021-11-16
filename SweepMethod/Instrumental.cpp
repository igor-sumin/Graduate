#include "Instrumental.h"


void Instrumental::printVec(const vec& a, const str& name) {
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

double Instrumental::getR(const vec& x, const vec& b) {
	int n = x.size();
	double R = 0.;

	for (int i = 0; i < n; i++) {
		R = fmax(R, fabs(b[i] - x[i]));
	}

	return R;
}

double Instrumental::getZ(const vec& u, const vec& v) {
	int n = u.size();
	double Z = 0.;

	for (int i = 0; i < n; i++) {
		Z = fmax(Z, fabs(u[i] - v[i]));
	}

	return Z;
}