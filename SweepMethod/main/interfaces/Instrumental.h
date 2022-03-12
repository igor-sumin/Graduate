#pragma once

#include <omp.h>

#include <iostream>
#include <vector>
#include <tuple>
#include <random>
#include <functional>
#include <algorithm>
#include <numeric>
#include "Task.h"

#define print() printf("---\n");
#define printd(a) printf("%s = %f\n", #a, a);
#define printi(a) printf("%s = %d\n", #a, a);

#define EPS 0.0001;

using matr = std::vector<std::vector<double>>;
using vec = std::vector<double>;
using pairs = std::pair<double, double>;
using str = std::string;


class Instrumental : public Task {
protected:
    size_t N;
    double h;
    vec x;

	vec u, v;
    TASK task;

    virtual void getDataByTask() {
        switch (task) {
            case TASK::TASK_7:
                this->defineDataByTask7();
                return;

            case TASK::NON_TASK:
                this->defineDataByNonTask();
                return;

            default:
                throw std::invalid_argument(Task::getValue(task));
        }
    }

public:
	Instrumental() : Instrumental(5) {}

	explicit Instrumental(size_t n) : Instrumental(n, TASK::NON_TASK) {}

    explicit Instrumental(size_t n, TASK task);

    // Getting a grid with nodes
    vec getGridNodes() const;

	// Printing a vector @a with @name
	static void printVec(const vec& a, const str& name);

	// Printing a matrix @a with @name
	static void printMatr(const matr& a, const str& name);

	// Calculating the discrepancy
	double calcR(const vec& x, const vec& b) const;

	// Calculating of the error estimate of the scheme
	double calcZ() const;

    const vec &getU() const;

    matr createThirdDiagMatr(double diag = 1., double upDiag = 3., double downDiag = 2.,
                             double first = 0.5, double second = 0.5,
                             double preLast = 0.5, double last = 0.5);

	// Matrix-vector multiplication : @A x @b
	static vec calcMatrVecMult(const matr& A, const vec& b);

	// Getting protected fields
	std::tuple<size_t, double, vec, vec, vec, TASK> getAllFields() const;

    // Setting protected fields
    void setAllFields(size_t N, double h,
                      const vec& x, const vec& u, const vec& v,
                      const TASK& task);

    // Compare of two doubles
    static bool compareDouble(double a, double b);

    // Comparison of two matrices
    static bool compareMatr(const matr& a, const matr& b);

    // Comparison of two vectors
    static bool compareVec(const vec& a, const vec& b);
};