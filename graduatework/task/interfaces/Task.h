#pragma once

#include <utility>

#include "System.h"

class Task : public System {
protected:
    // НУ
    InitConditions cond;

    // нелинейная функция F(u)
    vec3<vec3d<double>> F;

    // заданные функции для НУ
    // phi(x, y) = {phi[k](x, y) : k = 1,3 }
    vec2d<double> phi;

    vec2d<double> setConst(const vec3<double>& c) {
        vec2d<double> res;
        loop3(k) {
            res[k].assign(grid[0], vec(grid[1], c[k]));
        }

        return res;
    }

    vec2d<double> setCommonCos(const vec3<double>& m, double alpha, double beta) {
        vec2d<double> res;

        loop3(k) {
            for (size_t i = 0; i < grid[0]; i++) {
                for (size_t j = 0; j < grid[1]; j++) {
                    // res[k]_ij = M[k] * cos(pi * a * x_i) * cos(pi * b * y_j)
                    res[k][i][j] = m[k] * std::cos(pi * alpha * nodes[0][i]) * std::cos(pi * beta * nodes[1][j]);
                }
            }
        }

        return res;
    }

public:
    Task() = default;

    Task(InitConditions cond_, const Area& area_, const Grid& grid_, const Parameters& params_)
            : System(area_, grid_, params_), cond(std::move(cond_))
    {
        phi = definePhi();
        this->defineInitCond();
    }

    void defineInitCond() {
        // u[k](x, y, 0) = phi[k](x, y)
    }

    vec2d<double> definePhi() {
        auto setCos = [&](double alpha, double beta) {
            vec2d<double> C = setConst(cond.getC());
            vec2d<double> M = setCommonCos(cond.getM(), alpha, beta);

            // C + M * cosX * cosY
            std::transform(C.begin(), C.end(), M.begin(), M.begin(), std::plus<>());

            return C;
        };

        // define type of init cond
        switch(cond.getType()) {
            case Type::Const:
                phi = setConst(cond.getC());
                break;

            case Type::CosX:
                phi = setCos(cond.getAlpha(), 0.);
                break;

            case Type::CosY:
                phi = setCos(0., cond.getBeta());
                break;

            case Type::CosXCosY:
                phi = setCos(cond.getAlpha(), cond.getBeta());
                break;

            default:
                throw std::runtime_error("alarm4");
        }

        // phi[1](xi, yj) =
    }
};
