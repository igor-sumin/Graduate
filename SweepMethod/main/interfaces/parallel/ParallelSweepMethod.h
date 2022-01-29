#pragma once

#include <utility>

#include "main/interfaces/parallel/ParallelInstrumental.h"
#include "main/interfaces/AbstractSweepMethod.h"
#include "main/interfaces/serial/SerialSweepMethod.h"


class ParallelSweepMethod : public ParallelInstrumental, public AbstractSweepMethod {
protected:
	matr A;
	vec b, y;

	// Getting Y2 from vectors b, Y1
	vec getY2(vec& Y1);

	vec getX2(vec& Y2);

	vec result(const vec& X1, const vec& X2);

    // preprocessing Upper Left Corner of the matrix R
    void preULR(matr& R);

    // preprocessing Lower Right Corner of the matrix R
    void preLRR(matr& R);

public:
    ParallelSweepMethod() = default;

    ParallelSweepMethod(size_t n, size_t threadNum) : ParallelInstrumental(n, threadNum) {
        this->A = createThirdDiagMatrI();
        this->b = createVecN();
        this->y.assign(N, 0.);
    }

	ParallelSweepMethod(matr A_, vec b_) : A(std::move(A_)), b(std::move(b_)) {
        this->prepareData();
    }

    // Getting protected fields
    std::tuple<size_t, size_t, size_t, size_t, matr, vec, vec> getAllFields() const;

    // Setting protected fields
    void setAllFields(size_t N, size_t threadNum, size_t blockSize, size_t classicSize, matr A, vec b, vec y);

    void transformation();

    std::pair<matr, vec> collectInterferElem();

    vec collectPartY(const matr& R, const vec& partB);

    vec run() override;
};