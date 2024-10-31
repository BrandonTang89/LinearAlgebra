#include "LUFactorisation.h"
#include "Matrix.h"
#include "helpers.h"

#include <iostream>
bool LUFactorisation(const Matrix& A, Matrix& L, Matrix& U) {
    // Returns true if and only if the factorisation succeeded
    // Will POSSIBLY fail if the matrix is singular
    int n = A.H;
    if (A.W != n) {
        throw std::invalid_argument("Matrix must be square");
    }
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

                using namespace std;
                cout << U[{j, j}] << endl;
                if (std::abs(U[{j, j}]) < EPS) {
                    return false; // Matrix is singular
                }
                L[{i, j}] /= U[{j, j}];
                U[{i, j}] = 0.0;
            }
        }
    }
    return true; 
}
