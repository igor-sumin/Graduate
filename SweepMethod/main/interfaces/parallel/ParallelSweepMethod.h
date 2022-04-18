#pragma once

#include <utility>

#include "main/interfaces/parallel/ParallelInstrumental.h"
#include "main/interfaces/AbstractSweepMethod.h"
#include "main/interfaces/serial/SerialSweepMethod.h"


class ParallelSweepMethod : public ParallelInstrumental, public AbstractSweepMethod {
protected:

    // preprocessing Upper Left Corner of the matrix R
    void preULR(matr& R);

    // preprocessing Lower Right Corner of the matrix R
    void preLRR(matr& R);

public:
    ParallelSweepMethod() = default;

    ParallelSweepMethod(size_t n, size_t threadNum) : ParallelSweepMethod(n, threadNum, TASK::NON_TASK) {}

    ParallelSweepMethod(size_t n, size_t threadNum, TASK task) : ParallelInstrumental(n, threadNum, task) {}

	ParallelSweepMethod(matr A, vec b, vec y) : ParallelInstrumental(std::move(A), std::move(b), std::move(y)) {}

    /*
     * 1.
     * Reduction of the matrix to columnar views
    */
    void transformation();

    /*
     *  2.
     *  Collect the equations in the matrix- R and the right part - partB
     *  for the vector of unknown - partY, interfere with parallel execution
    */
    std::pair<matr, vec> collectInterferElem();

    /*
     * 3.
     * Solve of the system of linear algebraic equations
     * R * partY = partB by the classical sweep method
    */
    vec collectPartY(const matr& R, const vec& partB);

    /*
     * 4.
     * Find the elements of the vector of unknowns - y,
     * that do not interfere with parallel execution
    */
    void collectNotInterferElem();

    /*
     * 5.
     * Combine interfering elements - partY with part of non-interfering - y
    */
    vec collectFullY(const vec& partY);

    /*
     * Parallel implementation of sweep method
     *
     * N         - dimension
     * threadNum - number of compute nodes
     * blockSize - dimension for one computing node
     * interSize - dimension of interfering elements
     * A         - tridiagonal matrix (NxN)
     * b         - the right part of the system of linear algebraic equations
     * y         - vector of unknowns
    */
    vec run() override;
};