#include <gtest/gtest.h>

#include "../Matrix.h"
#include "../Vec.h"
#include "../ConjugateGradient.h"
#include <gtest/gtest.h>

#include "../Matrix.h"
#include "../Vec.h"
#include "../ConjugateGradient.h"

TEST(ConjugateGradientTest, Solve2x2SymmetricPositiveDefinite) {
    // Define a 2x2 symmetric positive definite matrix
    Matrix A(2, 2);
    A[{0, 0}] = 2.0;
    A[{0, 1}] = 1.0;
    A[{1, 0}] = 1.0;
    A[{1, 1}] = 3.0;

    // Define a 2x1 vector
    Vec b(2);
    b[0] = 1.0;
    b[1] = 2.0;

    // Solve the system Ax = b using Conjugate Gradient
    Vec x = conjugateGradient(A, b, 1000, 1e-16);
    
    Vec res = A * x - b;
    EXPECT_LT(res.norm(), 1e-9);
}
TEST(ConjugateGradientTest, Solve5x5SymmetricPositiveDefinite) {
    // Define a 5x5 lower triangular matrix L
    Matrix L(5, 5);
    L[{0, 0}] = 2.0;
    L[{1, 0}] = 3.0; L[{1, 1}] = 1.0;
    L[{2, 0}] = 1.0; L[{2, 1}] = 2.0; L[{2, 2}] = 1.0;
    L[{3, 0}] = 4.0; L[{3, 1}] = 1.0; L[{3, 2}] = 3.0; L[{3, 3}] = 1.0;
    L[{4, 0}] = 1.0; L[{4, 1}] = 2.0; L[{4, 2}] = 1.0; L[{4, 3}] = 2.0; L[{4, 4}] = 1.0;

    // Compute A = L * L^T to get a symmetric positive definite matrix
    Matrix A = L * L.transpose();

    // Define a 5x1 vector
    Vec b(5);
    b[0] = 1.0;
    b[1] = 2.0;
    b[2] = 3.0;
    b[3] = 4.0;
    b[4] = 5.0;

    // Solve the system Ax = b using Conjugate Gradient
    Vec x = conjugateGradient(A, b, 1000, 1e-16);
    
    Vec res = A * x - b;
    EXPECT_LT(res.norm(), 1e-9);
}

