#pragma once

#include <algo/interfaces/Instrumental.h>

#include <utility>
#include <ostream>

class Parameters {
public:
    // alpha1, alpha2
    pairs alpha;
    // beta1, beta2
    pairs beta;
    // gamma1, gamma2
    pairs gamma;

    Parameters() = default;

    Parameters(pairs alpha, pairs beta, pairs gamma) : alpha(std::move(alpha)), beta(std::move(beta)), gamma(std::move(gamma)) {}

    friend std::ostream &operator<<(std::ostream &os, const Parameters &parameters) {
        return os << "alpha: (" << parameters.alpha.first << ", " << parameters.alpha.second << "), "
                  << "beta: (" << parameters.beta.first << ", " << parameters.beta.second << "), "
                  << "gamma: (" << parameters.gamma.first << ", " << parameters.gamma.second << ")\n";
    }

    friend std::istream& operator>>(std::istream& in, Parameters& parameters) {
        return in >> parameters.alpha.first >> parameters.alpha.second
                  >> parameters.beta.first >> parameters.beta.second
                  >> parameters.gamma.first >> parameters.gamma.second;
    }
};