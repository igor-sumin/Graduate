#include "SerialInstrumental.h"


void SerialInstrumental::prepareData() {
}

bool SerialInstrumental::checkData() const {
	return true;
}

std::tuple<size_t, size_t, double, double, double, double, vec> SerialInstrumental::getAllFields() const {
	return std::make_tuple(N, node, h, A, C, B, u);
}


vec SerialInstrumental::getGridNodes() {
	vec x(node);

	for (int i = 0; i < node; i++) {
		x[i] = (double)i * h;
	}

	return x;
}