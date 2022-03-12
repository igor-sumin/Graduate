#include "main/interfaces/serial/SerialInstrumental.h"

void SerialInstrumental::setAllFields(const vec& A_, const vec& C_, const vec& B_,
                                      const vec& Phi_, pairs kappa_, pairs mu_, pairs gamma_) {
    this->A = A_;
    this->C = C_;
    this->B = B_;
    this->Phi = Phi_;
    this->kappa = kappa_;
    this->mu = mu_;
    this->gamma = gamma_;
}

std::tuple<vec, vec, vec, vec, pairs, pairs, pairs> SerialInstrumental::getAllFields() const {
    return std::make_tuple(A, C, B,
                           Phi,
                           kappa, mu, gamma);
}

void SerialInstrumental::prepareData() {
    if (!checkData()) {
        throw std::invalid_argument("n = " + std::to_string(N));
    }
}

bool SerialInstrumental::checkData() const {
    if (N < 0) {
        std::cerr << "Dimension (N = " << N << ") can't be negative.\n"
                << "Enter a positive dimension number (N)\n\n";

        return false;
    }

    return true;
}

void SerialInstrumental::defineDataByTask7() {
    double total = 12. / (h * h);

    A.assign(N - 1, total);
    C.assign(N - 1, 2. * total + 5.);
    B.assign(N - 1, total);

    mu = std::make_pair(10., 100.);
    kappa = std::make_pair(0., 0.);
    gamma = std::make_pair(1., 1.);

    Phi.assign(N - 1, 0.);
    for (size_t i = 0; i < N - 1; i++) {
        Phi[i] = 450. * x[i + 1] * x[i + 1] - 2110.;
    }
}

void SerialInstrumental::defineDataByNonTask() {
    A.assign(N - 1, 3.);
    C.assign(N - 1, 3.5);
    B.assign(N - 1, 3.);
    Phi.assign(N - 1, 1.);

    mu = std::make_pair(1., 1.);
    kappa = std::make_pair(0., 0.);
    gamma = std::make_pair(1., 1.);
}