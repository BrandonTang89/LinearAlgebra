#ifndef GAUSSIAN_ELIMINATION_H
#define GAUSSIAN_ELIMINATION_H

#include "Matrix.h"
#include "Vec.h"

double gaussianElimination(Matrix& A, Vec& ans, double& determinant);  // solve with A being the augmented matrix
double gaussianElimination(const Matrix& A, const Vec& b, Vec& ans);   // solve with A being the coefficient matrix
double determinant(const Matrix& A);

#endif  // GAUSSIAN_ELIMINATION_H