#include <cassert>
#include <cmath>
#include <iostream>

#include "Matrix.h"
#include "Vec.h"
Matrix getQ(const Matrix& A, int m, int n) {  // must have 1 <= m < M and 0 <= n
    const int M = A.H;
    const int N = A.W;
    assert(m < M && n < N);
    assert(m >= 1 && n >= 0);
    double denom = sqrt(A.data[m - 1][n] * A.data[m - 1][n] + A.data[m][n] * A.data[m][n]);
    double c = A.data[m - 1][n] / denom;
    double s = -A.data[m][n] / denom;
    Matrix Q = identity(M);
    Q.data[m - 1][m - 1] = c;
    Q.data[m][m] = c;
    Q.data[m - 1][m] = -s;
    Q.data[m][m - 1] = s;
    return Q;
}
std::pair<Matrix, Matrix> getQR(const Matrix& A) {  // not efficiently
    Matrix R = Matrix(A);
    int M = A.H;
    int N = A.W;
    Matrix Q = identity(M);
    for (int n = 0; n < N; n++) {
        for (int m = M - 1; m > n; m--) {
            Matrix Qmn = getQ(R, m, n);
            R = Qmn * R;
            Q = Q * transpose(Qmn);
        }
    }
    return {Q, R};
}
void givens(Matrix& A, Vec& B) {  // overwrites A with R and B with Qt R
    // Implementation of the algorithm in the lecture notes
    int M = A.H;
    int N = A.W;
    assert(B.N == M);
    for (int n = 0; n < N; n++) {
        for (int m = M - 1; m > n; m--) {
            double denom = sqrt(A.data[m - 1][n] * A.data[m - 1][n] + A.data[m][n] * A.data[m][n]);
            double c = A.data[m - 1][n] / denom;
            double s = -A.data[m][n] / denom;
            for (int k = n; k < N; k++) {
                double temp = A.data[m - 1][k];
                A.data[m - 1][k] = c * temp - s * A.data[m][k];
                A.data[m][k] = s * temp + c * A.data[m][k];
            }
            double temp = B[m - 1];
            B[m - 1] = c * temp - s * B[m];
            B[m] = s * temp + c * B[m];
        }
    }
}
// int main() {
//     std::cout << std::boolalpha;
//     int M = 5;
//     Matrix A(M, 2);
//     for (int i = 0; i < M; i++) {
//         A.data[i][0] = 1;
//         A.data[i][1] = i + 1;
//     }
//     std::cout << "A: " << std::endl;
//     std::cout << A << std::endl;
//     auto [Q, R] = getQR(A);
//     std::cout << "Q: " << Q << std::endl;
//     std::cout << "Q is orthogonal: " << isOrthogonal(Q) << std::endl;
//     std::cout << Q * transpose(Q) << std::endl;
//     std::cout << "R is upper triangular: " << isUpperTriangular(R) << std::endl;
//     std::cout << R << std::endl;
//     Matrix reconstructedA = Q * R;
//     std::cout << "Reconstructed A: " << std::endl;
//     std::cout << reconstructedA << std::endl;
//     std::cout << "A == reconstructed A: " << (A == reconstructedA) << std::endl;
//     std::cout << std::endl;
//     // Efficient Givens Rotation
//     Vec B(M, 11);
//     B[0] = 1;
//     std::cout << "B: " << B << std::endl;
//     std::cout << "Qt B: " << transpose(Q) * B << std::endl;

//     givens(A, B);
//     std::cout << "supposed R " << A << std::endl;
//     std::cout << "supposed QtB " << B << std::endl;
//     return 0;
// }