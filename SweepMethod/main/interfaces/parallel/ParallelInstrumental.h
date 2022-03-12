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
        this->getDataByTask();
    }

    ParallelInstrumental(matr A_, vec b_, vec y_) : A(std::move(A_)), b(std::move(b_)), y(std::move(y_)) {
        this->prepareData();
    }

    void setParallelOptions() const;

    // Getting protected fields
    std::tuple<size_t, size_t, size_t, size_t, matr, vec, vec> getAllFields() const;

    // Setting protected fields
    void setAllFields(size_t N, size_t threadNum, size_t blockSize, size_t classicSize, const matr& A_, const vec& b_, const vec& y_);

    // Preparing user data for parallel computing
    void prepareData(size_t threadNum);

    void prepareData() override;

    void defineDataByTask7() override;

    void defineDataByNonTask() override;

    // Checking for multiplicity of @N and @THREAD_NUM
    bool checkData() const override;

    // Creating a vector from 0 to @N with dimension @N
    vec createVecN();

    // Creating a vector right part for task 7
    vec createVecByTask7();
};