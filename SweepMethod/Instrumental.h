#pragma once

#include <iostream>
#include <omp.h>
#include <vector>
#include <tuple>
#include <random>
#include <algorithm>
#include <numeric>

using matr = std::vector<std::vector<double>>;
using vec = std::vector<double>;
using pairs = std::pair<double, double>;

class Instrumental {
private:
	bool isPrime(int num) const;
	
	// Preparing user data for parallel computing
	std::tuple<int, int, int, int> prepareData();
	
	// Ñhecking for multiplicity of @N and @THREADNUM
	bool checkData();
	
	// Finding all the divisors of @num
	vec findDivisors(int num);

protected:
	int N, THREAD_NUM, BLOCK_SIZE, CLASSIC_SIZE;

	Instrumental() : N(8), THREAD_NUM(4), BLOCK_SIZE(2), CLASSIC_SIZE(6) {}

	// Ñreating a vector from 0 to @N with dimension @N
	vec createVecN();

	// Creating a vector with random values from 0 to 100 with dimension @N
	vec createVecRand();

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

	/*
	 * Creating a tridiagonal matrix with dimension @N x @N
	 *
	 * side lower diagonal = @a
	 * side upper diagonal = @b
	 * main diagonal	   = @c
	*/
	matr createMatr(double a, double b, double c);

	// Matrix-vector multiplication : @A x @b
	vec calcMatrVecMult(const matr& A, const vec& b);

public:
	// Printing a vector @a with @name
	void printVec(const vec& a, const std::string& name);

	// Printing a matrix @a with @name
	void printMatr(const matr& a, const std::string& name);
};