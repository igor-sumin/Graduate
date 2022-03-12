#pragma once

#include "main/interfaces/Instrumental.h"

#include <utility>


class SerialInstrumental : public Instrumental {
protected:
    vec A, C, B;

    vec Phi;
    pairs kappa, mu, gamma;

public:
    SerialInstrumental() = default;

    explicit SerialInstrumental(size_t n) : SerialInstrumental(n, TASK::NON_TASK) {}

    // TODO: -1 ?
    SerialInstrumental(size_t n, TASK task) : Instrumental(n - 1, task) {
        this->prepareData();
        this->getDataByTask();
    }

    SerialInstrumental(vec a, vec c, vec b, vec phi, pairs kappa_, pairs mu_, pairs gamma_) :
        SerialInstrumental(std::move(a), std::move(c), std::move(b), std::move(phi), kappa_, mu_, gamma_, TASK::NON_TASK)
    {}

    SerialInstrumental(vec a, vec c, vec b, vec phi, pairs kappa_, pairs mu_, pairs gamma_, TASK task_) :
    // TODO: +1 ?
        Instrumental(a.size() + 1, task_),
            A(std::move(a)), C(std::move(c)), B(std::move(b)),
            Phi(std::move(phi)),
            kappa(std::move(kappa_)), mu(std::move(mu_)), gamma(std::move(gamma_))
    {
        this->prepareData();
    }

    // Preparing user data for serial computing
    void prepareData() override;

    bool checkData() const override;

    void defineDataByTask7() override;

    void defineDataByNonTask() override;

    // Setting protected fields
    void setAllFields(const vec& A, const vec& C, const vec& B,
                      const vec& Phi_, pairs kappa_, pairs mu_, pairs gamma_);

    // Getting protected fields
    std::tuple<vec, vec, vec, vec, pairs, pairs, pairs> getAllFields() const;
};