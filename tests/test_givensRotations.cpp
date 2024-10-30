#include <gtest/gtest.h>
#include "../GivensRotations.h"
#include "../Matrix.h"
#include "../Vec.h"

TEST(GivensRotationsTest, QRDecomposition) {
    int M = 5;
    Matrix A(M, 2);
    for (int i = 0; i < M; i++) {
        A.data[i][0] = 1;
        A.data[i][1] = i + 1;
    }

    auto [Q, R] = getQR(A);

    EXPECT_TRUE(isOrthogonal(Q));
    EXPECT_TRUE(isUpperTriangular(R));

    Matrix reconstructedA = Q * R;
    EXPECT_EQ(A, reconstructedA);
}

TEST(GivensRotationsTest, EfficientGivensRotation) {
    int M = 5;
    Matrix A(M, 2);
    for (int i = 0; i < M; i++) {
        A.data[i][0] = 1;
        A.data[i][1] = i + 1;
    }

    Vec B(M, 11);
    B[0] = 1;

    givens(A, B);

    // Assuming we have functions to check if A is upper triangular and B is transformed correctly
    EXPECT_TRUE(isUpperTriangular(A));
    // Should add checks that B = Qt * B
}