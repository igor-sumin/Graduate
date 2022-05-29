#pragma once

#include <vector>
#include <utility>

#include "task/interfaces/constants/AppConstants.h"
#include "algo/interfaces/serial/SerialSweepMethod.h"

// x, y, t
template<class T>
using vec3 = std::vector<T>;
template<class T>
using matr2d = std::vector<std::vector<T>>;
template<class T>
using vec2d = vec3<matr2d<T>>;
template<class T>
using matr3d = std::vector<std::vector<std::vector<T>>>;
template<class T>
using vec3d = vec3<matr3d<T>>;

#define loop3(x) for (size_t x = 0; (x) < 3; (x)++)
#define loop(x, n) for (size_t x = 0; (x) < (n); (x)++)

#include "task/interfaces/tools/Parameters.h"
#include "task/interfaces/tools/Area.h"
#include "task/interfaces/tools/Grid.h"
#include "InitConditions.h"
#include "fstreams/interfaces/FWork.h"

class System {
private:
    void assignAllU() {
        matr2d<double> null(grid[0], vec(grid[1], 0.));

        uLow.assign(3, null);
        uMid.assign(3, null);
        uTop.assign(3, null);
    }

public:
    // обл опред функции
    vec3<pairs> area;

    // сетка
    // n1, n2, m
    vec3<int> grid;

    // параметры для функции F
    // alpha12, beta12, omega12, gamma12
    vec3<pairs> params;

    // шаг сетки
    // h1, h2, tau
    vec3<double> step;

    // узлы сетки
    // x_i, y_j, t_s
    vec3<vec> nodes;

    // точное решение ду
    // u(x, y, t) = { u[k](x, y, t) : k = 1,3 }
    // vec3<matr3d<double>> u;

    // используются следующие слои:
    // -> текущий       - low
    // -> промежуточный - mid
    // -> новый         - top
    vec2d<double> uLow;
    vec2d<double> uMid;
    vec2d<double> uTop;

    // подвижность
    vec3<double> lambda;

//    // коэф-ты 1 и 2 слоев (A1, A2)
//    vec3<pairs> A, B, C;
//
//    // правая часть слоев
//    vec2d<pairs> Phi;

    // коэф-ты 1 и 2 слоев (A1, A2)
    vec3<double> A, B, C;

    // правая часть слоев
    vec2d<double> Phi;

    // нелинейная функция F(u)
    vec2d<double> F;

    // точное решение разностной схемы
    // v(x, y, t) = { v[k](x, y, t) : k = 1,3 }
    // vec3<matr3d<double>> v;

    System() = default;

    System(const Area& area_, const Grid& grid_, Parameters params_,
           const vec3<double>& lambda_)
            : area(area_.getData()), grid(grid_.getData()), params(params_.getData()),
              lambda(lambda_)
    {
        this->assignAllU();
        this->defineStep();
        this->defineNodes();

        Phi.assign(3, matr2d<double>(grid[0], vec(grid[1], 0.)));
    }

    void defineStep() {
        step.resize(3);
        loop3(i) {
            step[i] = (area[i].second - area[i].first) / (grid[i] - 1);
        }
    }

    void defineNodes() {
        nodes.resize(3);
        loop3(k) {
            nodes[k].assign(grid[k], 0.);

            loop(ijs, grid[k]) {
                nodes[k][ijs] = area[k].first + (double)ijs * step[k];
            }
        }
    }

    void defineF(size_t phase) {
        F = (phase == 0)
                ? this->declareF(uLow)
                : this->declareF(uMid);
    }

    vec2d<double> declareF(const vec2d<double>& u) {
        vec2d<double> res(3, matr2d<double>(grid[0], vec(grid[1], 0.)));

        for (int i = 0; i < grid[0]; i++) {
            for (int j = 0; j < grid[1]; j++) {
                res[0][i][j] = u[0][i][j] * (params[0].first - u[0][i][j] - params[2].first * u[1][i][j] - params[1].first * u[2][i][j]);
                res[1][i][j] = u[1][i][j] * (params[0].second - params[2].second * u[0][i][j] - u[1][i][j] - params[1].second * u[2][i][j]);
                res[2][i][j] = u[2][i][j] * (-1. + params[3].first * u[0][i][j] + params[3].second * u[1][i][j]);
            }
        }

        return res;
    }

    void defineLayersParams(size_t phase) {
        A.assign(3, 0.);
        B.assign(3, 0.);
        C.assign(3, 0.);

        // для 1 и 2 слоев
        loop3(k) {
            // tau * lambda_k / h12^2
            auto common = [&](double h12) {
                return step[2] * lambda[k] / std::pow(h12, 2);
            };

            // -0.5 * tau * lambda_k / h12^2
            auto formulaAB = [&](double h12) {
                return -0.5 * common(h12);
            };

            // -(1 + tau * lambda_k / h12^2)
            auto formulaC = [&](double h12) {
                return (-1) * (1 + common(h12));
            };

            // для правой части
            auto layer = [&](size_t i0, size_t j0, size_t i1, size_t j1, size_t i2, size_t j2,
                             double h12, const vec2d<double>& u) {
                return (-1) * (u[k][i0][j0] * (1 - common(h12)) + 0.5 * common(h12) * (u[k][i1][j1] + u[k][i2][j2]) + F[k][i0][j0]);
            };


            A[k] = B[k] = (phase == 0)
                            ? formulaAB(step[0])
                            : formulaAB(step[1]);

            C[k] = (phase == 0)
                        ? -formulaC(step[0])
                        : -formulaC(step[1]);

            if (phase == 0) {
                for (size_t i = 1; i < grid[0] - 1; i++) {
                    for (size_t j = 1; j < grid[1] - 1; j++) {
                        Phi[k][i][j] = layer(i, j, i, j - 1, i, j + 1, step[1], uLow);
                    }
                }

            } else if (phase == 1) {
                for (size_t i = 1; i < grid[0] - 1; i++) {
                    for (size_t j = 1; j < grid[1] - 1; j++) {
                        Phi[k][i][j] = layer(i, j, i - 1, j, i + 1, j, step[0], uMid);
                    }
                }

            } else {
                throw std::runtime_error(AppConstansts::ALARM_LAYERS_PARAMS);
            }
        }
    }

    void executeSerialSweepLayersPhase(size_t phase, vec2d<double>& uPhase) {
        size_t n1 = (phase == 0) ? grid[0] : grid[1];
        size_t n2 = (phase == 0) ? grid[1] : grid[0];

        pairs kappa = std::make_pair(1., 1.);
        pairs mu = std::make_pair(0., 0.);
        pairs gamma = std::make_pair(1., 1.);

        loop3(k) {
            vec a(n1 - 2, A[k]);
            vec c(n1 - 2, C[k]);
            vec b(n1 - 2, B[k]);
            vec phi(n1 - 2, 0.);

            for (size_t ind1 = 0; ind1 < n2; ind1++) {
                for (size_t ind2 = 1; ind2 < n1 - 1; ind2++) {
                    // -Phi[k]_1js ... -Phi[k]_n1-1js
                    phi[ind2 - 1] = Phi[k][ind2][ind1];
                }

                SerialSweepMethod ssm(a, c, b, phi, kappa, mu, gamma);
                vec u = ssm.run();

                for (size_t ind2 = 0; ind2 < n1; ind2++) {
                    uPhase[k][ind2][ind1] = u[ind2];
                }
            }
        }
    }

    void executeSerialSweepLayers(size_t phase) {
        if (phase == 0) {
            this->executeSerialSweepLayersPhase(phase, uMid);

        } else if (phase == 1) {
            this->executeSerialSweepLayersPhase(phase, uTop);

        } else {
            throw std::runtime_error(AppConstansts::ALARM_EXECUTE_LAYERS);
        }
    }
};