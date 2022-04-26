#include <test/ParallelAlgorithmComponentTest.h>
#include <test/SerialAlgorithmComponentTest.h>
#include <test/common/InstrumentalComponentTest.h>

#include "task/impl/Task.h"

int main(int argc, char** argv) {
//	ParallelAlgorithmComponentTest p;
//    p.execute();
//
//    SerialAlgorithmComponentTest s;
//    s.execute();
//
//    InstrumentalComponentTest i;
//    i.execute();

    Task task;
    int n1, n2, m;
    Parameters params;

    // задаем размер сетки
    std::cout << "Enter (n1, n2, m):\n";
    std::cin >> n1 >> n2 >> m;

    // задаем параметры для задачи
    std::cout << "Enter (alpha1, alpha2), (beta1, beta2), (gamma1, gamma2):\n";
    std::cin >> params;

    // задаем

	return 0;
}