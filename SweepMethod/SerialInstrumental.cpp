#include "SerialInstrumental.h"


void SerialInstrumental::prepareData() {
	this->getGridNodes();
}

bool SerialInstrumental::checkData() const {
	return true;
}

std::tuple<vec, double, double, double, double> SerialInstrumental::getAllFields() {
	this->prepareData();
	this->checkData();

	return std::make_tuple(x, h, A, C, B);
}


vec SerialInstrumental::getGridNodes() {
	for (int i = 0; i < node; i++) {
		x[i] = (double)i * h;
	}

	return x;
}

matr SerialInstrumental::createMatr() {
	matr res(node, vec(node));

	print(A);
	print(B);
	print(C);

	#pragma omp parallel for if (node > 500)
	for (int i = 1; i < node - 1; i++) {
		for (int j = 0; j < node; j++) {
			res[i][i] = C;
			res[i][i - 1] = A;
			res[i - 1][i] = B;
		}
	}

	res[0][0] = 1.; res[node - 1][node - 1] = 1.;
	res[0][1] = 0.; res[node - 1][node - 2] = 0.;

	return res;
}