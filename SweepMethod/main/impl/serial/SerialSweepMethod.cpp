#include "main/interfaces/serial/SerialSweepMethod.h"


vec SerialSweepMethod::run() {
    vec alpha(N), beta(N);
    vec y(N);

    // forward motion
    alpha[0] = kappa.first / gamma.first;
    beta[0] = mu.first / gamma.first;
    for (size_t i = 0; i < N - 1; i++) {
        alpha[i + 1] = B[i] / (C[i] - A[i] * alpha[i]);
        beta[i + 1] = (Phi[i] + A[i] * beta[i]) / (C[i] - A[i] * alpha[i]);
    }

    // reverse motion
    y[N - 1] = (-kappa.second * beta[N - 1] - mu.second) / (kappa.second * alpha[N - 1] - gamma.second);
    for (int i = (int)N - 2; i >= 0; i--) {
        y[i] = alpha[i] * y[i + 1] + beta[i];
    }

    return y;
}