#pragma once

#include <utility>

#include "../../algo/interfaces/Instrumental.h"
#include "Parameters.h"
#include "Area.h"
#include "Grid.h"

// x, y, t
template<class T>
using vec3 = std::vector<T>;
template<class T>
using matr3d = std::vector<std::vector<std::vector<T>>>;
template<class T>
using vec3d = vec3<matr3d<T>>;

class System {
public:
    // обл опред функции
    Area area;

    // сетка
    Grid grid;

    // параметры для функции F
    Parameters params;

    // шаг сетки
    vec3<double> step;

    // узлы сетки
    vec3<vec> nodes;

    // точное решение ду
    // u(x, y, t) = { u[k](x, y, t) : k = 1,3 }
    vec3<matr3d<double>> u;

    // точное решение разностной схемы
    // v(x, y, t) = { v[k](x, y, t) : k = 1,3 }
    vec3<matr3d<double>> v;

    // заданные функции для НУ
    // phi(x, y) = {phi[k](x, y) : k = 1,3 }
    vec3<std::function<double(double, double)>> phi;

};