#include "main/interfaces/serial/SerialSweepMethod.h"

#include <utility>


vec SerialSweepMethod::run() {
    vec alpha(N), beta(N);
    vec y(node);

    // forward motion
    alpha[0] = kappa.first / gamma.first;
    beta[0] = mu.first / gamma.first;
    for (size_t i = 0; i < N - 1; i++) {
         alpha[i + 1] = B[i] / (C[i] - A[i] * alpha[i]);
         beta[i + 1] = (Phi[i] + A[i] * beta[i]) / (C[i] - A[i] * alpha[i]);
    }

    // reverse motion
    y[node - 1] = (-kappa.second * beta[N - 1] - mu.second) / (kappa.second * alpha[N - 1] - gamma.second);
    for (int i = (int)N - 1; i >= 0; i--) {
        y[i] = alpha[i] * y[i + 1] + beta[i];
    }

    return y;
}

std::tuple<vec, vec, size_t, size_t, double, vec, vec, vec, vec, vec, pairs, pairs, pairs> SerialSweepMethod::getAllFields() {
    return std::make_tuple(v, u, N, node, h, x, A, C, B, Phi, kappa, mu, gamma);
}

void SerialSweepMethod::setAllFields(vec v, vec u, size_t N, size_t node, double h, vec x, vec A, vec C, vec B, vec Phi, pairs kappa, pairs mu, pairs gamma) {
    this->v = std::move(v);
    this->u = std::move(u);
    this->N = N;
    this->node = node;
    this->h = h;
    this->x = std::move(x);
    this->A = std::move(A);
    this->C = std::move(C);
    this->B = std::move(B);
    this->Phi = std::move(Phi);
    this->kappa = kappa;
    this->mu = mu;
    this->gamma = gamma;
}