#include "SerialSweepMethod.h"


/*
 * Последовательная реализация метода прогонки
 *
 * b - диагональ, лежащая над главной
 * c - главная диагональ матрицы А
 * а - дигональ, лежащая под главной
 * f - правая часть
 * x - вектор решений
 * kappa, mu, gamma - коэффициенты
 */
vec SerialSweepMethod::run() {
	vec alpha(n - 1), beta(n - 1);
	vec x(n);

	// forward motion
	alpha[0] = kappa.first / gamma.first;
	beta[0] = mu.first / gamma.first;;
	for (size_t i = 0; i < n - 2; i++) {
		alpha[i + 1] = b[i] / (c[i] - a[i] * alpha[i]);
		beta[i + 1] = (f[i] + a[i] * beta[i]) / (c[i] - a[i] * alpha[i]);
	}

	// reverse motion
	x[n - 1] = (-kappa.second * beta[n - 2] - mu.second) / (kappa.second * alpha[n - 2] - gamma.second);
	for (size_t i = n - 2; i >= 0; i--) {
		x[i] = alpha[i] * x[i + 1] + beta[i];
	}

	return x;
}