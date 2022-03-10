#pragma once

#include "main/interfaces/Instrumental.h"


class ParallelInstrumental : public Instrumental {
private:
    // is the @num prime
    static bool isPrime(int num);

    // Finding all the divisors of @num
    static vec findDivisors(int num);

protected:
    matr A;
    vec b, y;

    size_t threadNum, blockSize, interSize;

public:
    ParallelInstrumental() = default;

    ParallelInstrumental(size_t n, size_t tN) : ParallelInstrumental(n, tN, TASK::NON_TASK) {}

    ParallelInstrumental(size_t n, size_t tN, TASK task) : Instrumental(n, task) {
        this->prepareData(tN);
        this->defineDataByTask();
    }

    ParallelInstrumental(matr A_, vec b_) : A(std::move(A_)), b(std::move(b_)) {
        this->prepareData();
        this->defineDataByTask();
    }

    void setParallelOptions() const;

    // Preparing user data for parallel computing
    void prepareData(size_t threadNum);

    void prepareData() override;

    void defineDataByTask7() override;

    void defineDataByNonTask() override;

    // Checking for multiplicity of @N and @THREAD_NUM
    bool checkData() const override;

    matr createThirdDiagMatr(double diag, double upDiag, double downDiag,
                             double first, double second,
                             double preLast, double last);

    // Creating a vector from 0 to @N with dimension @N
    vec createVecN();

    // Creating a vector right part for task 7
    vec createVecByTask7();

    // Creating a vector of result for task 7
    vec createResByTask7();
};