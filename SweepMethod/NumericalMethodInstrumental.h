#pragma once

// Получить сетку с узлами
vec getX(int n) {
	vec x(n);

	double h = 1 / static_cast<double>(n);
	for (int i = 0; i < n; i++) {
		x[i] = i * h;
	}

	return x;
}

// Подсчет невязки
double getR(const vec& x, const vec& b) {
	int n = x.size();
	double R = 0.;

	for (int i = 0; i < n; i++) {
		R = fmax(R, fabs(b[i] - x[i]));
	}

	return R;
}

// Подсчет оценки погрешности схемы
double getZ(const vec& u, const vec& v) {
	int n = u.size();
	double Z = 0.;

	for (int i = 0; i < n; i++) {
		Z = fmax(Z, fabs(u[i] - v[i]));
	}

	return Z;
}