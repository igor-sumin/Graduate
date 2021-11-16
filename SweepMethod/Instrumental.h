#pragma once

#include <omp.h>

#include <iostream>
#include <vector>
#include <tuple>
#include <random>
#include <functional>
#include <algorithm>
#include <numeric>

#define print(a) printf("%s = %f\n", #a, a)

using matr = std::vector<std::vector<double>>;
using vec = std::vector<double>;
using pairs = std::pair<double, double>;
using str = std::string;


class Instrumental {
protected:
	size_t N;

public:
	Instrumental() : N(5) {}

	Instrumental(size_t n) : N(n) {}

	// Preparing user data for parallel computing
	virtual void prepareData() = 0;

	// Checking for multiplicity of @N and @THREADNUM
	virtual bool checkData() const = 0;

	// Printing a vector @a with @name
	static void printVec(const vec& a, const str& name);

	// Printing a matrix @a with @name
	static void printMatr(const matr& a, const str& name);

	// Calculating the discrepancy
	static double getR(const vec& x, const vec& b);

	// Calculating of the error estimate of the scheme
	static double getZ(const vec& u, const vec& v);
};