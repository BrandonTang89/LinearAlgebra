#include <iostream>
#include <random>

#include "ConjugateGradient.h"
#include "Matrix.h"
#include "Vec.h"

int main() {

    std::random_device rd;   // obtain a random number from hardware
    std::mt19937 gen(rd());  // seed the generator
                                                    
    // Create a uniform real distribution between 0 and 1
    std::uniform_real_distribution<> dis(0.0, 1.0);

    int N = 10;
    Matrix L = Matrix(N, N);  // a lower triangular matrix
    for (int i = 0; i < N; i++) {
        L.data[i][i] = static_cast<double>(i + 1);
        for (int j = 0; j < i; j++) {
            L.data[i][j] = dis(gen);
        }
    }

    Matrix A = L * transpose(L);  // A is a symmetric positive definite matrix
    Vec b = Vec(N, 0);
    for (int i = 0; i < N; i++) {
        b[i] = dis(gen);
    }

    Vec x = conjugateGradient(A, b, 1000, 1e-16);
    std::cout << "Solution to Ax = b using Conjugate Gradient: " << x << std::endl;
    std::cout << "Residual: " << (A * x - b).norm() << std::endl;

    // std::cout << "A = " << A << std::endl;
    Vec x2 = conjugateGradientPreconditionedDiagonal(A, b, 1000, 1e-16);
    std::cout << "Solution to Ax = b using Preconditioned Conjugate Gradient: " << x2 << std::endl;
    std::cout << "Residual: " << (A * x2 - b).norm() << std::endl;
    return 0;
}