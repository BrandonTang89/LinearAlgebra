#include "Vec.h"
#include "Matrix.h"
#include <cmath>

// Copy assignment operator
Vec& Vec::operator=(const Vec& v) {
    if (this != &v) {
        N = v.N;
        data = v.data;
    }
    return *this;
}

// Stream insertion operator
std::ostream& operator<<(std::ostream& os, const Vec& v) {
    os << "Vector of size " << v.N << ":\n";
    for (int i = 0; i < v.N; ++i) {
        os << v.data[i] << " ";
    }
    os << std::endl;
    return os;
}

// Transpose function
Matrix Vec::transpose() const {
    Matrix result(N, 1);
    for (int i = 0; i < N; ++i) {
        result.data[i][0] = data[i];
    }
    return result;
}

// Norm function
double Vec::norm() {
    double sum = 0.0;
    for (double val : data) {
        sum += val * val;
    }
    return std::sqrt(sum);
}

// Vector addition
Vec Vec::operator+(const Vec& other) const {
    if (N != other.N) throw std::invalid_argument("Vector dimensions do not match for addition.");
    Vec result(N);
    for (int i = 0; i < N; ++i) {
        result[i] = data[i] + other[i];
    }
    return result;
}

// Vector subtraction
Vec Vec::operator-(const Vec& other) const {
    if (N != other.N) throw std::invalid_argument("Vector dimensions do not match for subtraction.");
    Vec result(N);
    for (int i = 0; i < N; ++i) {
        result[i] = data[i] - other[i];
    }
    return result;
}

// Scalar multiplication
Vec operator*(double c, const Vec& v) {
    Vec result(v.N);
    for (int i = 0; i < v.N; ++i) {
        result[i] = c * v[i];
    }
    return result;
}

// Dot Product
double Vec::dot(const Vec& other) const {
    if (N != other.N) throw std::invalid_argument("Vector dimensions do not match for dot product.");
    double result = 0.0;
    for (int i = 0; i < N; ++i) {
        result += data[i] * other[i];
    }
    return result;
}