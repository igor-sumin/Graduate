#pragma once

#include "ParallelInstrumental.h"
#include "AbstractSweepMethod.h"
#include "SerialSweepMethod.h"

class ParallelSweepMethod : public ParallelInstrumental, public AbstractSweepMethod {
private:
	matr A;
	vec b;

	std::pair<matr, vec> collectInterferElem();

	vec ClassicSweepMethod(const matr& R, const vec& Y1);

	// Получаем Y2 из векторов b, Y1
	vec getY2(vec& Y1);

	vec getX2(vec& Y2);

	vec result(const vec& X1, const vec& X2);

	void transformMatrA();

public:
	ParallelSweepMethod() : A(createThirdDiagMatrI()), b(createVecN()) {}
	ParallelSweepMethod(matr A_, vec b_) : A(A_), b(b_) {}

	vec run() override;
};