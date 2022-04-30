#pragma once

#include <set>
#include <cmath>

const double pi = std::acos(-1);

enum class Type {
    NotSelect, Const, CosX, CosY, CosXCosY

    //
};

class InitConditions {
private:
    Type type;

    vec3<double> C;
    vec3<double> M;
    double alpha, beta;

public:

    InitConditions() = default;

    InitConditions(Type type) : type(type) {}

    Type getType() const {
        return type;
    }

    vec3<double> getC() const {
        return C;
    }

    vec3<double> getM() const {
        return M;
    }

    double getAlpha() const {
        return alpha;
    }

    double getBeta() const {
        return beta;
    }

    friend std::ostream &operator<<(std::ostream &os, const InitConditions &conditions) {
        str c = "(C1, C2, C3)";
        str m = "(M1, M2, M3)";
        str alp = "alpha (= 1,2,3..)";
        str bet = "beta (= 1,2,3..)";

        switch(conditions.getType()) {
            case Type::Const:
                return os << "Enter" << c << ":\n";

            case Type::CosX:
                return os << "Enter" << c << ", " << m << ", " << alp << ":\n";

            case Type::CosY:
                return os << "Enter" << c << ", " << m << ", " << bet << ":\n";

            case Type::CosXCosY:
                return os << "Enter" << c << ", " << m << ", " << alp << ", " << bet << ":\n";

            case Type::NotSelect: default:
                throw std::runtime_error("alarm1");
        }
    }

    friend std::istream &operator>>(std::istream &in, InitConditions &conditions) {
        auto getConst = [&](vec3<double>& Const) {
            loop3(k) {
                in >> Const[k];
            }
        };

        switch(conditions.getType()) {
            case Type::Const:
                getConst(conditions.C);
                break;

            case Type::CosX:
                getConst(conditions.C);
                getConst(conditions.M);
                in >> conditions.alpha;
                break;

            case Type::CosY:
                getConst(conditions.C);
                getConst(conditions.M);
                in >> conditions.beta;
                break;

            case Type::CosXCosY:
                getConst(conditions.C);
                getConst(conditions.M);
                in >> conditions.alpha >> conditions.beta;
                break;

            case Type::NotSelect: default:
                throw std::runtime_error("alarm3");
        }

        return in;
    }
};

str getValue(const Type& type)  {
    switch(type) {
        case Type::NotSelect:
            return "not selected yet";
        case Type::Const:
            return "(1) const";
        case Type::CosX:
            return "(2) const + cos(x)";
        case Type::CosY:
            return "(3) const + cos(y)";
        case Type::CosXCosY:
            return "(4) const + cos(x)cos(y)";
    }

    ASSERT(false);
    return {};
}

std::ostream &operator<<(std::ostream &os, const Type &type) {
    if (type == Type::NotSelect) {
        return os << "Choose one of four possible types of initial conditions:\n"
                  << "(1) const\n" << "(2) const + cos(x)\n" << "(3) const + cos(y)\n"
                  << "(4) const + cos(x)cos(y)\n"
                  << "Enter 1 or 2 or 3 or 4:\n";
    }

    return os << "your type is - " << getValue(type);
}

std::istream &operator>>(std::istream &in, Type &type) {
    std::set<int> input = {1, 2, 3, 4};
    int output;

    do {
        in >> output;

        switch(output) {
            case 1:
                type = Type::Const;
                break;
            case 2:
                type = Type::CosX;
                break;
            case 3:
                type = Type::CosY;
                break;
            case 4:
                type = Type::CosXCosY;
                break;
            default:
                std::cerr << "Enter 1 or 2 or 3 or 4:\n";
        }
    } while (input.find(output) == input.end());

    return in;
}