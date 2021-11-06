#include <iostream>
#include <fstream>
#include <omp.h>
#include <vector>
#include <tuple>
#include <string>
#include <random>
#include <algorithm>
#include <numeric>
#include <functional>
#include <cmath>
#include "stdio.h"

using matr = std::vector<std::vector<double>>;
using vec = std::vector<double>;
using pairs = std::pair<double, double>;

#include "Instrumental.h"
#include "NumericalMethodInstrumental.h"

#include "SerialSweepMethod.h"
#include "ParallelSweepMethod.h"

#include "UnitTests.h"

int main(int argc, char** argv) {
	setlocale(LC_ALL, "rus");
	callAllTests();
	
	return 0;
}