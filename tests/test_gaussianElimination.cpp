#include <gtest/gtest.h>

#include "../Matrix.h"
#include "../Vec.h"
#include "../GaussianElimination.h"

TEST(GaussianEliminationTest, NoSolution) {
    Matrix A(2, 2);
    A[{0, 0}] = 1.0;
    A[{0, 1}] = 1.0;
    A[{1, 0}] = 1.0;
    A[{1, 1}] = 1.0;

    Vec b(2);
    b[0] = 1.0;
    b[1] = 2.0;

    Vec ans(2);
    int result = gaussianElimination(A, b, ans);

    EXPECT_EQ(result, 0);
}

TEST(GaussianEliminationTest, OneSolution) {
    Matrix A(2, 2);
    A[{0, 0}] = 1.0;
    A[{0, 1}] = 1.0;
    A[{1, 0}] = 1.0;
    A[{1, 1}] = 2.0;

    Vec b(2);
    b[0] = 1.0;
    b[1] = 2.0;

    Vec ans(2);
    int result = gaussianElimination(A, b, ans);

    EXPECT_EQ(result, 1);
    EXPECT_DOUBLE_EQ(ans[0], 0.0);
    EXPECT_DOUBLE_EQ(ans[1], 1.0);
}

TEST(GaussianEliminationTest, InfiniteSolutions) {
    Matrix A(2, 2);
    A[{0, 0}] = 1.0;
    A[{0, 1}] = 1.0;
    A[{1, 0}] = 2.0;
    A[{1, 1}] = 2.0;

    Vec b(2);
    b[0] = 1.0;
    b[1] = 2.0;

    Vec ans(2);
    int result = gaussianElimination(A, b, ans);

    EXPECT_EQ(result, 2);
}

TEST(GaussianEliminationTest, Determinant0) {
    Matrix A(2, 2);
    A[{0, 0}] = 1.0;
    A[{0, 1}] = 1.0;
    A[{1, 0}] = 1.0;
    A[{1, 1}] = 1.0;

    double det = determinant(A);

    EXPECT_DOUBLE_EQ(det, 0.0);
}

TEST(GaussianEliminationTest, DeterminantPositive) {
    Matrix A(2, 2);
    A[{0, 0}] = 1.0;
    A[{0, 1}] = 1.0;
    A[{1, 0}] = 1.0;
    A[{1, 1}] = 2.0;

    double det = determinant(A);

    EXPECT_DOUBLE_EQ(det, 1.0);
}
TEST(GaussianEliminationTest, NegativeDeterminant) {
    Matrix A(2, 2);
    A[{0, 0}] = 1.0;
    A[{0, 1}] = 2.0;
    A[{1, 0}] = 3.0;
    A[{1, 1}] = 4.0;

    double det = determinant(A);

    EXPECT_DOUBLE_EQ(det, -2.0);
}
