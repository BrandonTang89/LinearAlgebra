#include <gtest/gtest.h>
#include "../LUFactorisation.h"
#include "../Matrix.h"

TEST(LUFactorisationTest, RegularMatrix) {
    Matrix A(3, 3);
    A[{0, 0}] = 2.0; A[{0, 1}] = -1.0; A[{0, 2}] = -2.0;
    A[{1, 0}] = -4.0; A[{1, 1}] = 6.0; A[{1, 2}] = 3.0;
    A[{2, 0}] = -4.0; A[{2, 1}] = -2.0; A[{2, 2}] = 8.0;

    Matrix L, U;
    bool isSingular = LUFactorisation(A, L, U);

    EXPECT_FALSE(isSingular);

    EXPECT_TRUE(isLowerTriangular(L));
    EXPECT_TRUE(isUpperTriangular(U));
    EXPECT_EQ(L * U, A);
}

TEST(LUFactorisationTest, SingularMatrix) {
    Matrix A(3, 3);
    A[{0, 0}] = 1.0; A[{0, 1}] = 2.0; A[{0, 2}] = 3.0;
    A[{1, 0}] = 4.0; A[{1, 1}] = 5.0; A[{1, 2}] = 6.0;
    A[{2, 0}] = 7.0; A[{2, 1}] = 8.0; A[{2, 2}] = 9.0;

    Matrix L, U;
    bool isSingular = LUFactorisation(A, L, U);

    EXPECT_TRUE(isSingular);
}
