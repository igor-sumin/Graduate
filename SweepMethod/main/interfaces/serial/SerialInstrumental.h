#pragma once

#include "main/interfaces/Instrumental.h"


class SerialInstrumental : public Instrumental {
protected:
    vec x;
    double h;
    double A, C, B;

public:
    SerialInstrumental() : SerialInstrumental(5) {}

    explicit SerialInstrumental(size_t n) :
            Instrumental(n),
            x(node), h(1 / static_cast<double>(N)), A(-1), C(-1), B(-1) {}

    // Preparing user data for serial computing
    void prepareData() override;

    // Checking for ...
    bool checkData() const override;

    // Getting a grid with nodes
    vec getGridNodes();

    // Getting protected fields
    std::tuple<vec, double, double, double, double> getAllFields();

    /*
     * Creating a tridiagonal matrix with dimension @N x @N
     *
     * side lower diagonal = @a
     * side upper diagonal = @b
     * main diagonal	   = @c
    */
    matr createMatr();
};