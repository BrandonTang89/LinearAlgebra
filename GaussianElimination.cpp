#include "GaussianElimination.h"

#include "Matrix.h"
#include "Vec.h"
#include "helpers.h"

double gaussianElimination(Matrix& A, Vec& ans, double& determinant) {
    // returns 0 if no solution, 1 if unique solution, 2 if infinite solutions
    // also returns the determinant of the matrix
    int n = A.H;
    int m = A.W - 1;

    std::vector<int> where(m, -1);
    determinant = 1;
    for (int col = 0, row = 0; col < m && row < n; ++col) {
        int sel = row;
        for (int i = row; i < n; ++i) {
            if (std::abs(A[{i, col}]) > std::abs(A[{sel, col}])) {
                sel = i;
            }
        }
        if (std::abs(A[{sel, col}]) < EPS) {
            determinant = 0;
            continue;
        }

        if (sel != row) {
            determinant *= -1;
            for (int i = col; i <= m; ++i) {
                std::swap(A[{sel, i}], A[{row, i}]);
            }
        }

        determinant *= A[{row, col}];
        double scalar = 1.0 / A[{row, col}];  // normalise the row
        for (int i = col; i <= m; ++i) {
            A[{row, i}] *= scalar;
        }

        where[col] = row;
        for (int i = 0; i < n; ++i) {
            if (i != row) {
                double c = A[{i, col}];
                for (int j = col; j <= m; ++j) {
                    A[{i, j}] -= A[{row, j}] * c;
                }
            }
        }
        ++row;
    }

    ans = Vec(m);
    for (int i = 0; i < m; ++i) {
        if (where[i] != -1) {
            ans[i] = A[{where[i], m}] / A[{where[i], i}];
        }
    }
    for (int i = 0; i < n; ++i) {
        double sum = 0;
        for (int j = 0; j < m; ++j) {
            sum += ans[j] * A[{i, j}];
        }
        if (std::abs(sum - A[{i, m}]) > EPS) {
            return 0;  // No solution
        }
    }

    for (int i = 0; i < m; ++i) {
        if (where[i] == -1) {
            return 2;  // Infinite solutions
        }
    }
    return 1;  // Unique solution
}

double gaussianElimination(const Matrix& A, const Vec& b, Vec& ans) {
    Matrix augmented(A.H, A.W + 1);
    for (int i = 0; i < A.H; ++i) {
        for (int j = 0; j < A.W; ++j) {
            augmented[{i, j}] = A[{i, j}];
        }
        augmented[{i, A.W}] = b[i];
    }

    double determinant;
    return gaussianElimination(augmented, ans, determinant);
}

double determinant(const Matrix& A) {
    Matrix augmented(A.H, A.W + 1);
    for (int i = 0; i < A.H; ++i) {
        for (int j = 0; j < A.W; ++j) {
            augmented[{i, j}] = A[{i, j}];
        }
    }

    Vec ans(A.H);
    double determinant;
    gaussianElimination(augmented, ans, determinant);
    return determinant;
}