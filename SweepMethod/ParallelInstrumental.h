#pragma once

#include "Instrumental.h"


class ParallelInstrumental : public Instrumental {
private:
	bool isPrime(int num) const;

	// Finding all the divisors of @num
	vec findDivisors(int num) const;

protected:
	int THREAD_NUM, BLOCK_SIZE, CLASSIC_SIZE;

public:
	ParallelInstrumental();

	ParallelInstrumental(size_t n, int threadNum, int blockSize, int classicSize)
		: Instrumental(n), 
		  THREAD_NUM(threadNum), BLOCK_SIZE(blockSize), CLASSIC_SIZE(classicSize) {}

	// Preparing user data for parallel computing
	void prepareData() override;

	// Checking for multiplicity of @N and @THREADNUM
	bool checkData() const override;

	// Getting protected fields
	std::tuple<size_t, int, int, int> getAllFields() const;

	/*
	 * Creating a tridiagonal matrix with dimension @N x @N
	 *
	 * side lower diagonal = @a
	 * side upper diagonal = @b
	 * main diagonal	   = @c
	*/
	matr createMatr(double a, double b, double c);

	/*
	 * Creating a tridiagonal matrix with dimension @N x @N
	 *
	 * side lower diagonal = 1
	 * side upper diagonal = 2
	 * main diagonal	   = 3
	*/
	matr createThirdDiagMatrI();

	// Creating a matrix with random values from 0 to 100 with dimension @N x @N
	matr createThirdDiagMatrRand();

	// Matrix-vector multiplication : @A x @b
	vec calcMatrVecMult(const matr& A, const vec& b);

	// Creating a vector from 0 to @N with dimension @N
	vec createVecN();

	// Creating a vector with random values from 0 to 100 with dimension @N
	vec createVecRand();
};