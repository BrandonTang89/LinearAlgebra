#ifndef CONJUGATEGRADIENT_H
#define CONJUGATEGRADIENT_H

#include <functional>
#include "Vec.h"

Vec conjugateGradient(const Matrix& A, const Vec& b, int maxIter, double tol);
Vec conjugateGradientPreconditioned(const Matrix& A, const Vec& b, int maxIter, double tol, std::function<Vec(Vec)> multMInv);
Vec conjugateGradientPreconditionedDiagonal(const Matrix& A, const Vec& b, int maxIter, double tol); 

#endif // CONJUGATEGRADIENT_H