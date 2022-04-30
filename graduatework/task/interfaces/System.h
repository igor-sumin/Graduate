#pragma once

#include <utility>

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

class System {
public:
    // обл опред функции
    vec3<pairs> area;

    // сетка
    // n1, n2, m
    vec3<int> grid;

    // параметры для функции F
    // alpha12, beta12, gamma12
    vec3<pairs> params;

    // шаг сетки
    // h1, h2, tau
    vec3<double> step;

    // узлы сетки
    // x_i, y_j, t_s
    vec3<vec> nodes;

    // точное решение ду
    // u(x, y, t) = { u[k](x, y, t) : k = 1,3 }
    vec3<matr3d<double>> u;

    // точное решение разностной схемы
    // v(x, y, t) = { v[k](x, y, t) : k = 1,3 }
    vec3<matr3d<double>> v;

    System() = default;

    System(const Area& area_, const Grid& grid_, Parameters params_)
        : area(area_.getData()), grid(grid_.getData()), params(params_.getData())
    {
        this->defineStep();
        this->defineNodes();
    }

    void defineStep() {
        step.assign(3, 0.);

        loop3(i) {
            step[i] = (area[i].second - area[i].first) / grid[i];
        }
    }

    void defineNodes() {
        nodes.assign(3, vec());

        loop3(k) {
            nodes[k].assign(grid[k], 0.);

            loop(ijs, grid[k]) {
                nodes[k][ijs] = area[k].first + (double)ijs * step[k];
            }
        }
    }
};