#pragma once

#include <ostream>

class Grid {
public:
    int n1, n2, m;

    Grid() = default;

    friend ostream &operator<<(ostream &os, const Grid &grid) {
        return os << "n1: " << grid.n1 << " n2: " << grid.n2 << " m: " << grid.m;
    }

    friend std::istream& operator>>(std::istream& in, Grid& grid) {
        return in >> grid.n1 >> grid.n2 >> grid.m;
    }
};