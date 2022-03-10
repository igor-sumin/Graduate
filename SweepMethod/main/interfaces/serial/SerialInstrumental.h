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

    SerialInstrumental(size_t n, TASK task) : Instrumental(n, task) {
        this->prepareData();
        this->defineDataByTask();
    }

    // Preparing user data for serial computing
    void prepareData() override;

    bool checkData() const override;

    void defineDataByTask7() override;

    void defineDataByNonTask() override;
};