#include <gtest/gtest.h>

#include "../Matrix.h"
#include "../Vec.h"

// Test the Vec norm function
TEST(VecTest, Norm) {
    Vec vec(3);
    vec[0] = 1.0;
    vec[1] = 2.0;
    vec[2] = 2.0;
    EXPECT_EQ(vec.norm(), 3.0);
}

// Test the Vec transpose function
TEST(VecTest, Transpose) {
    Vec vec(3);
    vec[0] = 1.0;
    vec[1] = 2.0;
    vec[2] = 3.0;
    Matrix transposed = vec.transpose();
    Matrix expected(1, 3);
    expected[{0, 0}] = 1.0;
    expected[{0, 1}] = 2.0;
    expected[{0, 2}] = 3.0;
    EXPECT_EQ(transposed.H, expected.H);
    EXPECT_EQ(transposed.W, expected.W);

    EXPECT_EQ(transposed, expected);
}
// Test the Vec addition function
TEST(VecTest, Addition) {
    Vec vec1(3);
    vec1[0] = 1.0;
    vec1[1] = 2.0;
    vec1[2] = 3.0;

    Vec vec2(3);
    vec2[0] = 4.0;
    vec2[1] = 5.0;
    vec2[2] = 6.0;

    Vec result = vec1 + vec2;
    EXPECT_EQ(result[0], 5.0);
    EXPECT_EQ(result[1], 7.0);
    EXPECT_EQ(result[2], 9.0);
}

// Test the Vec subtraction function
TEST(VecTest, Subtraction) {
    Vec vec1(3);
    vec1[0] = 4.0;
    vec1[1] = 5.0;
    vec1[2] = 6.0;

    Vec vec2(3);
    vec2[0] = 1.0;
    vec2[1] = 2.0;
    vec2[2] = 3.0;

    Vec result = vec1 - vec2;
    EXPECT_EQ(result[0], 3.0);
    EXPECT_EQ(result[1], 3.0);
    EXPECT_EQ(result[2], 3.0);
}

// Test the Vec scalar multiplication function
TEST(VecTest, ScalarMultiplication) {
    Vec vec(3);
    vec[0] = 1.0;
    vec[1] = 2.0;
    vec[2] = 3.0;

    Vec result = 2.0 * vec;
    EXPECT_EQ(result[0], 2.0);
    EXPECT_EQ(result[1], 4.0);
    EXPECT_EQ(result[2], 6.0);
}
