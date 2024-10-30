#ifndef GIVENS_ROTATIONS_H
#define GIVENS_ROTATIONS_H

#include "Matrix.h"
#include "Vec.h"

Matrix getQ(const Matrix& A, int m, int n);
std::pair<Matrix, Matrix> getQR(const Matrix& A);
void givens(Matrix& A, Vec& B);

#endif  // GIVENS_ROTATIONS_H