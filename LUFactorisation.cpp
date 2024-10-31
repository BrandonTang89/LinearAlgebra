#include "LUFactorisation.h"
#include "Matrix.h"
#include "helpers.h"

bool LUFactorisation(const Matrix& A, Matrix& L, Matrix& U) {
    int n = A.H;
    L = Matrix(n, n);
    U = Matrix(n, n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i <= j) {
                U[{i, j}] = A[{i, j}];
                for (int k = 0; k < i; ++k) {
                    U[{i, j}] -= L[{i, k}] * U[{k, j}];
                }
                if (i == j) {
                    L[{i, j}] = 1.0;
                } else {
                    L[{i, j}] = 0.0;
                }
            } else {
                L[{i, j}] = A[{i, j}];
                for (int k = 0; k < j; ++k) {
                    L[{i, j}] -= L[{i, k}] * U[{k, j}];
                }
                if (std::abs(U[{j, j}]) < EPS) {
                    return false; // Matrix is singular
                }
                L[{i, j}] /= U[{j, j}];
                U[{i, j}] = 0.0;
            }
        }
    }
    return true; // Matrix is not singular
}
