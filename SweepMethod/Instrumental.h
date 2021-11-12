#pragma once

#include <omp.h>

#include <iostream>
#include <vector>
#include <tuple>
#include <random>
#include <functional>
#include <algorithm>
#include <numeric>

using matr = std::vector<std::vector<double>>;
using vec = std::vector<double>;
using pairs = std::pair<double, double>;
using str = std::string;


class Instrumental {
private:
	bool isPrime(int num) const;
	
	// Finding all the divisors of @num
	vec findDivisors(int num);

protected:
	int N, THREAD_NUM, BLOCK_SIZE, CLASSIC_SIZE;

public:
	// Instrumental() : N(8), THREAD_NUM(4), BLOCK_SIZE(2), CLASSIC_SIZE(6) {}
	Instrumental();

	Instrumental(int n, int threadNum, int blockSize, int classicSize) 
		: N(n), THREAD_NUM(threadNum), BLOCK_SIZE(blockSize), CLASSIC_SIZE(classicSize) {}

	// Preparing user data for parallel computing
	void prepareData();

	// Checking for multiplicity of @N and @THREADNUM
	bool checkData();

	// Printing a vector @a with @name
	static void printVec(const vec& a, const str& name);

	// Printing a matrix @a with @name
	static void printMatr(const matr& a, const str& name);

	// Getting a grid with nodes
	vec getX();

	// Calculating the discrepancy
	double getR(const vec & x, const vec & b);

	// Calculating of the error estimate of the scheme
	double getZ(const vec & u, const vec & v);

	// Getting protected fields
	std::tuple<int, int, int, int> getFields() const;

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