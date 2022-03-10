#include "main/interfaces/serial/SerialInstrumental.h"


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

    for (size_t i = 0; i < node - 2; i++) {
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