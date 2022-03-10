#pragma once

#include <utility>

#include "main/interfaces/serial/SerialInstrumental.h"
#include "main/interfaces/AbstractSweepMethod.h"


class SerialSweepMethod : public SerialInstrumental, public AbstractSweepMethod {
public:
    SerialSweepMethod() = default;

    explicit SerialSweepMethod(size_t n) : SerialSweepMethod(n, TASK::NON_TASK){}

    SerialSweepMethod(size_t n, TASK task) : SerialInstrumental(n, task) {}


    SerialSweepMethod(vec a, vec c, vec b, vec phi, pairs kappa_, pairs mu_, pairs gamma_) :
        SerialInstrumental(std::move(a), std::move(c), std::move(b), std::move(phi), kappa_, mu_, gamma_) {}

    // Getting protected fields
    std::tuple<vec, vec, size_t, size_t, double, vec,
               vec, vec, vec,
               vec, pairs, pairs, pairs> getFields() const;

    // Setting protected fields
    void setAllFields(const vec& v, const vec& u, size_t N, size_t node, double h, const vec& x,
                      const vec& A, const vec& C, const vec& B,
                      const vec& Phi_, pairs kappa_, pairs mu_, pairs gamma_);

    vec method();

    /*
     * Sequential implementation of sweep method
     *
     * a   - is the diagonal lying under the main
     * c   - the main diagonal of the matrix A
     * b   - the diagonal lying above the main
     * phi - right side
     * x   - vector of solutions
     * kappa, mu, gamma - coefficients
    */
    vec run() override;
};