#pragma once

#include <utility>

#include "main/interfaces/parallel/ParallelInstrumental.h"
#include "main/interfaces/AbstractSweepMethod.h"
#include "main/interfaces/serial/SerialSweepMethod.h"


class ParallelSweepMethod : public ParallelInstrumental, public AbstractSweepMethod {
protected:
	matr A;
	vec b, y;

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
    void setAllFields(size_t N, size_t threadNum, size_t blockSize, size_t classicSize, const matr& A_, const vec& b_, const vec& y_);

    void transformation();

    std::pair<matr, vec> collectInterferElem();

    vec collectPartY(const matr& R, const vec& partB);

    void collectNotInterferElem();

    void collectFullY(const vec& partY);

    /*
     * ...
    */
    vec run() override;
};