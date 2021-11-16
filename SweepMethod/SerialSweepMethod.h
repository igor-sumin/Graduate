#pragma once

#include "SerialInstrumental.h"
#include "AbstractSweepMethod.h"


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
						: n(c_.size() - 2), 
						  a(a_), b(b_), f(f_), 
						  kappa(kappa_), mu(mu_), gamma(gamma_) {}

	vec run() override;
};
