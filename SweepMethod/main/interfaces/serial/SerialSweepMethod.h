#pragma once

#include <utility>

#include "main/interfaces/serial/SerialInstrumental.h"
#include "main/interfaces/AbstractSweepMethod.h"


class SerialSweepMethod : public SerialInstrumental, public AbstractSweepMethod {
private:
	int n;
	vec a, c, b, f;
	pairs kappa, mu, gamma;

public:
	SerialSweepMethod() : n(8),
						  a(10, 1), c(10, 2), b(10, 3), f(10, 4),
						  kappa(std::make_pair(0., 0.)), mu(std::make_pair(1., 1.)), gamma(std::make_pair(1., 1.)) {}

	SerialSweepMethod(vec a_, vec c_, vec b_, vec f_, pairs kappa_, pairs mu_, pairs gamma_)
						: n((int)c_.size() - 2),
						  a(std::move(a_)), c(std::move(c_)), b(std::move(b_)), f(std::move(f_)),
						  kappa(std::move(kappa_)), mu(std::move(mu_)), gamma(std::move(gamma_)) {}

    /*
     * Sequential implementation of sweep method
     *
     * b - the diagonal lying above the main
     * c - the main diagonal of the matrix A
     * a - is the diagonal lying under the main
     * f - right side
     * x - vector of solutions
     * kappa, mu, gamma - coefficients
    */
	vec run() override;
};
