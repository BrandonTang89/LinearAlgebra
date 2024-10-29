#include <gtest/gtest.h>
#include "../Matrix.h"

// Test the Matrix constructor
TEST(MatrixTest, Constructor) {
    Matrix mat(3, 3);
    EXPECT_EQ(mat.H, 3);
    EXPECT_EQ(mat.W, 3);
    EXPECT_EQ(mat.data.size(), 3);
    EXPECT_EQ(mat.data[0].size(), 3);
}

// Test the Matrix copy constructor
TEST(MatrixTest, CopyConstructor) {
    Matrix mat1(3, 3);
    mat1[{0, 0}] = 1.0;
    mat1[{1, 1}] = 2.0;
    mat1[{2, 2}] = 3.0;

    Matrix mat2 = mat1;
    EXPECT_EQ(mat2, mat1);
}


// Test the Matrix multiplication
TEST(MatrixTest, MultiplicationDirectComparison) {
    Matrix mat1(2, 3);
    mat1[{0, 0}] = 1.0; mat1[{0, 1}] = 2.0; mat1[{0, 2}] = 3.0;
    mat1[{1, 0}] = 4.0; mat1[{1, 1}] = 5.0; mat1[{1, 2}] = 6.0;

    Matrix mat2(3, 2);
    mat2[{0, 0}] = 7.0; mat2[{0, 1}] = 8.0;
    mat2[{1, 0}] = 9.0; mat2[{1, 1}] = 10.0;
    mat2[{2, 0}] = 11.0; mat2[{2, 1}] = 12.0;

    Matrix expected(2, 2);
    expected[{0, 0}] = 58.0; expected[{0, 1}] = 64.0;
    expected[{1, 0}] = 139.0; expected[{1, 1}] = 154.0;

    Matrix result = mat1 * mat2;

    EXPECT_EQ(result, expected);
}

// Test the Matrix addition
TEST(MatrixTest, AdditionDirectComparison) {
    Matrix mat1(2, 2);
    mat1[{0, 0}] = 1.0; mat1[{0, 1}] = 2.0;
    mat1[{1, 0}] = 3.0; mat1[{1, 1}] = 4.0;

    Matrix mat2(2, 2);
    mat2[{0, 0}] = 5.0; mat2[{0, 1}] = 6.0;
    mat2[{1, 0}] = 7.0; mat2[{1, 1}] = 8.0;

    Matrix expected(2, 2);
    expected[{0, 0}] = 6.0; expected[{0, 1}] = 8.0;
    expected[{1, 0}] = 10.0; expected[{1, 1}] = 12.0;

    Matrix result = mat1 + mat2;

    EXPECT_EQ(result, expected);
}

// Test the Matrix scaling
TEST(MatrixTest, ScalingDirectComparison) {
    Matrix mat(2, 2);
    mat[{0, 0}] = 1.0; mat[{0, 1}] = 2.0;
    mat[{1, 0}] = 3.0; mat[{1, 1}] = 4.0;

    Matrix expected(2, 2);
    expected[{0, 0}] = 2.0; expected[{0, 1}] = 4.0;
    expected[{1, 0}] = 6.0; expected[{1, 1}] = 8.0;

    Matrix result = 2.0 * mat;

    EXPECT_EQ(result, expected);
}


