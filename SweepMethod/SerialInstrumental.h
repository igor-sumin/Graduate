#pragma once

#include "Instrumental.h"


class SerialInstrumental : public Instrumental {
protected:
	size_t node;
	double h;
	double A, C, B;
	vec u;

public:
	SerialInstrumental() : SerialInstrumental(5) {}

	SerialInstrumental(size_t n) :
		Instrumental(n),
		node(N + 1), h(1 / static_cast<double>(N)), A(-1), C(-1), B(-1), u(node) {}

	// Preparing user data for parallel computing
	void prepareData() override;

	// Checking for multiplicity of @N and @THREADNUM
	bool checkData() const override;

	// Getting a grid with nodes
	vec getGridNodes();

	// Getting protected fields
	std::tuple<size_t, size_t, double, double, double, double, vec> getAllFields() const;
};