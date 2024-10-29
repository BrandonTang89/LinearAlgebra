#ifndef CONJUGATEGRADIENT_H
#define CONJUGATEGRADIENT_H

#include "Vec.h"

Vec conjugateGradient(const Matrix& A, const Vec& b, int maxIter, double tol);
#endif // CONJUGATEGRADIENT_H