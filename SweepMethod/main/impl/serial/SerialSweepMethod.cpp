#include "main/interfaces/serial/SerialSweepMethod.h"


vec SerialSweepMethod::run() {
	vec alpha(n - 1), beta(n - 1);
	vec x(n);

	// forward motion
	alpha[0] = kappa.first / gamma.first;
	beta[0] = mu.first / gamma.first;;
	for (int i = 0; i < n - 2; i++) {
		alpha[i + 1] = b[i] / (c[i] - a[i] * alpha[i]);
		beta[i + 1] = (f[i] + a[i] * beta[i]) / (c[i] - a[i] * alpha[i]);
	}

	// reverse motion
	x[n - 1] = (-kappa.second * beta[n - 2] - mu.second) / (kappa.second * alpha[n - 2] - gamma.second);
	for (int i = n - 2; i >= 0; i--) {
		x[i] = alpha[i] * x[i + 1] + beta[i];
	}

	return x;
}