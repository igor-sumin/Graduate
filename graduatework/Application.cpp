#include <task/interfaces/Task.h>
#include <test/ParallelAlgorithmComponentTest.h>

int main(int argc, char** argv) {
    Area area(std::make_pair(0, 1.6), std::make_pair(0, 1.6), std::make_pair(0, 1));
    Grid grid(16, 16, 100);
    Parameters params(std::make_pair(3, 1), std::make_pair(2, 4), std::make_pair(4, 1), std::make_pair(20, 1));
    Type type = Type::CosXCosY;
    // Type type = Type::Const;
    vec lambdas = {1., 1., 200.};

    // задаем размер сетки
    std::cout << "Enter (n1, n2, m):\n";

    // задаем параметры для задачи
    std::cout << "Enter (alpha1, alpha2), (beta1, beta2), (omega1, omega2), (gamma1, gamma2):\n";
    std::cout << "Enter (lambda1, lambda2, lambda3):\n";

    // задаем размер области
    std::cout << "Enter [a, b] x [c, d] x [0, T]:\n";

    // задаем тип НУ
    std::cout << type;

    // InitConditions cond(type, {0.05, 0.19, 0.65}, {0.01, -0.001, 0.002}, 1, 1);
     InitConditions cond(type, {0.04, 0.2, 0.64}, {0.01, -0.02, 0.03}, 1, 1);
//    InitConditions cond(type, {0.05, 0.21, 0.65}, {0.01, -0.02, 0.03}, 1, 1);
    // InitConditions cond(type, {0.04, 0.2, 0.64});

    // задаем НУ
    std::cout << cond;

    Task task(cond, area, grid, params, lambdas, true);
    task.execute();

    return 0;
}