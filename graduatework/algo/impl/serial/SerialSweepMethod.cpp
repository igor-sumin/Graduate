#include "algo/interfaces/serial/SerialSweepMethod.h"

#include <utility>


std::tuple<vec, vec, size_t, size_t, double, vec, vec, vec, vec, vec, pairs, pairs, pairs> SerialSweepMethod::getFields() const {
    return std::make_tuple(v, u,
                           N, node, h, x,
                           A, C, B,
                           Phi, kappa, mu, gamma);
}

void SerialSweepMethod::setAllFields(const vec& v, const vec& u, size_t N, size_t node, double h, const vec& x,
                                     const vec& A, const vec& C, const vec& B,
                                     const vec& Phi_, pairs kappa_, pairs mu_, pairs gamma_) {
    this->v = v;
    this->u = u;
    this->N = N;
    this->node = node;
    this->h = h;
    this->x = x;
    this->A = A;
    this->C = C;
    this->B = B;
    this->Phi = Phi_;
    this->kappa = kappa_;
    this->mu = mu_;
    this->gamma = gamma_;
}

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