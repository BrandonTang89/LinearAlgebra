#include <gtest/gtest.h>
#include "../LUFactorisation.h"
#include "../Matrix.h"
#include "../GaussianElimination.h"
#include <iostream>

TEST(LUFactorisationTest, RegularMatrix) {
    Matrix A(3, 3);
    A[{0, 0}] = 2.0; A[{0, 1}] = -1.0; A[{0, 2}] = -2.0;
    A[{1, 0}] = -4.0; A[{1, 1}] = 6.0; A[{1, 2}] = 3.0;
    A[{2, 0}] = -4.0; A[{2, 1}] = -2.0; A[{2, 2}] = 8.0;

    Matrix L(3, 3);
    Matrix U(3, 3);
    bool isInvertible = LUFactorisation(A, L, U);

    EXPECT_TRUE(isInvertible);
    EXPECT_TRUE(isLowerTriangular(L));
    EXPECT_TRUE(isUpperTriangular(U));
    EXPECT_EQ(L * U, A);
}
