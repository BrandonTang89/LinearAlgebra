#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>
#include <tuple>

struct Vec;

struct Matrix {
    int H;
    int W;
    std::vector<std::vector<double>> data;

    // Constructors
    explicit Matrix(const Vec& v);
    Matrix(int h, int w);
    Matrix(const Matrix& A);
    Matrix(Matrix&& A) noexcept;

    // Assignment Operators
    Matrix& operator=(const Matrix& A);
    Matrix& operator=(Matrix&& A) noexcept;

    // Matrix addition operator
    Matrix operator+(const Matrix& B) const;
    Matrix operator-(const Matrix& B) const;

    // Matrix multiplication operator
    Matrix operator*(const Matrix& B) const;

    // Matrix-Vector multiplication operator
    Vec operator*(const Vec v) const;

    // Matrix equality operator
    bool operator==(const Matrix& B) const;

    // Matrix Scaling
    friend Matrix operator*(double scalar, const Matrix& A);

    // Stream insertion operator
    friend std::ostream& operator<<(std::ostream& os, const Matrix& A);

    // Matrix Transpose
    Matrix transpose() const;

    // C++23 multi-subscripts operator[]
    double& operator[](std::tuple<int, int> rc){
        return data[std::get<0>(rc)][std::get<1>(rc)];
    }
    double operator[](std::tuple<int, int> rc) const {
        return data[std::get<0>(rc)][std::get<1>(rc)];
    }

};

// Transpose function
Matrix transpose(const Matrix& A);
Matrix identity(int n);

bool isOrthogonal(const Matrix& Q);
bool isUpperTriangular(const Matrix& A);
bool isLowerTriangular(const Matrix& A);

#endif // MATRIX_H
