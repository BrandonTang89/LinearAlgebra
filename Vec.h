#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <vector>

// Forward declaration of Matrix class
struct Matrix;

struct Vec {
    int N{};
    std::vector<double> data;
    Vec(int n) : N(n), data(n, 0) {}
    Vec(int n, double val) : N(n), data(n, val) {}
    Vec(const Vec& v) : N(v.N), data(v.data) {}
    Vec& operator=(const Vec& v);  // copy assignment operator

    // Random access operator
    double& operator[](int i) { return data[i]; }
    const double& operator[](int i) const { return data[i]; }

    friend std::ostream& operator<<(std::ostream& os, const Vec& v);

    // Vector addition
    Vec operator+(const Vec& other) const;
    Vec operator-(const Vec& other) const;

    // Scalar multiplication
    friend Vec operator*(double c, const Vec& v);

    Matrix transpose() const;
    double norm();
    double dot(const Vec& other) const;
};

#endif  // VECTOR_H