#pragma once

#include <algo/interfaces/Instrumental.h>

#include <utility>

// G U dG = [a, b] x [c, d] x [0, T]
class Area {
public:
    // [a, b]
    pairs x;
    // [c, d]
    pairs y;
    // [0, T]
    pairs t;

    Area() = default;

    Area(pairs x, pairs y, pairs t) : x(std::move(x)), y(std::move(y)), t(std::move(t)) {}

    friend std::ostream &operator<<(std::ostream &os, const Area &parameters) {
        return os << "x: (" << parameters.x.first << ", " << parameters.x.second << "), "
                  << "y: (" << parameters.y.first << ", " << parameters.y.second << "), "
                  << "t: (" << parameters.t.first << ", " << parameters.t.second << ")\n";
    }

    friend std::istream& operator>>(std::istream& in, Area& parameters) {
        return in >> parameters.x.first >> parameters.x.second
                  >> parameters.y.first >> parameters.y.second
                  >> parameters.t.first >> parameters.t.second;
    }
};