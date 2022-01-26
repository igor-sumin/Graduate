#pragma once

#include <utility>

#include "main/interfaces/parallel/ParallelInstrumental.h"
#include "main/interfaces/AbstractSweepMethod.h"
#include "main/interfaces/serial/SerialSweepMethod.h"


class ParallelSweepMethod final : public ParallelInstrumental, public AbstractSweepMethod {
private:
	matr A;
	vec b;

	std::pair<matr, vec> collectInterferElem();

	vec ClassicSweepMethod(const matr& R, const vec& Y1);

	// Getting Y2 from vectors b, Y1
	vec getY2(vec& Y1);

	vec getX2(vec& Y2);

	vec result(const vec& X1, const vec& X2);

public:
    ParallelSweepMethod(size_t n, size_t threadNum) : ParallelInstrumental(n, threadNum) {
        this->A = createThirdDiagMatrI();
        this->b = createVecN();
    }

	ParallelSweepMethod(matr A_, vec b_) : A(std::move(A_)), b(std::move(b_)) {
        this->prepareData();
    }

    // Getting protected fields
    std::tuple<size_t, size_t, size_t, size_t, matr, vec> getAllFields() const;

    // Setting protected fields
    void setAllFields(size_t N, size_t threadNum, size_t blockSize, size_t classicSize, matr A, vec b);

    void transformation();

    vec run() override;
};