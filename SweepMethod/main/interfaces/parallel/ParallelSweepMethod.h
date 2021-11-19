#pragma once

#include "main/interfaces/tools/Instrumental.h"
#include "main/interfaces/AbstractSweepMethod.h"
#include "main/interfaces/serial/SerialSweepMethod.h"

class ParallelSweepMethod : public Instrumental, public AbstractSweepMethod {
private:
	matr A;
	vec b;

	std::pair<matr, vec> collectInterferElem();

	vec ClassicSweepMethod(const matr& R, const vec& Y1);

	// Getting Y2 from vectors b, Y1
	vec getY2(vec& Y1);

	vec getX2(vec& Y2);

	vec result(const vec& X1, const vec& X2);

	void transformMatrA();

public:
	ParallelSweepMethod() : A(createThirdDiagMatrI()), b(createVecN()) {}
	ParallelSweepMethod(matr A_, vec b_) : A(A_), b(b_) {}

	vec run() override;
};