#include <matplot/matplot.h>

#include <iostream>
#include <random>

// #include "ConjugateGradient.h"
// #include "Matrix.h"
// #include "Vec.h"
#include "gillespie.h"

// int main() {

//     std::random_device rd;   // obtain a random number from hardware
//     std::mt19937 gen(rd());  // seed the generator

//     // Create a uniform real distribution between 0 and 1
//     std::uniform_real_distribution<> dis(0.0, 1.0);

//     int N = 10;
//     Matrix L = Matrix(N, N);  // a lower triangular matrix
//     for (int i = 0; i < N; i++) {
//         L.data[i][i] = static_cast<double>(i + 1);
//         for (int j = 0; j < i; j++) {
//             L.data[i][j] = dis(gen);
//         }
//     }

//     Matrix A = L * transpose(L);  // A is a symmetric positive definite matrix
//     Vec b = Vec(N, 0);
//     for (int i = 0; i < N; i++) {
//         b[i] = dis(gen);
//     }

//     Vec x = conjugateGradient(A, b, 1000, 1e-16);
//     std::cout << "Solution to Ax = b using Conjugate Gradient: " << x << std::endl;
//     std::cout << "Residual: " << (A * x - b).norm() << std::endl;

//     // std::cout << "A = " << A << std::endl;
//     Vec x2 = conjugateGradientPreconditionedDiagonal(A, b, 1000, 1e-16);
//     std::cout << "Solution to Ax = b using Preconditioned Conjugate Gradient: " << x2 << std::endl;
//     std::cout << "Residual: " << (A * x2 - b).norm() << std::endl;
//     return 0;
// }

int main() {
    using namespace matplot;

    double A0 = 1000;
    double k1 = 0.1;
    Reaction dup = {
        [=](std::vector<double> reactants, double volume) { return k1 * reactants[0]; },
        {1}};

    double k2 = 0.01;
    Reaction decay = {
        [=](std::vector<double> reactants, double volume) { return k2 * reactants[0]; },
        {-1}};

    std::vector<Reaction> reactions = {decay, dup};
    std::vector<double> reactants = {A0};
    auto [times, trajectory] = glllespie(50, reactants, reactions, 1.0);
    double nTimePoints = trajectory.size();
    std::vector<double> a_amt(nTimePoints);
    for (int i = 0; i < nTimePoints; i++) {
        a_amt[i] = trajectory[i][0];
    }

    auto compMean = [=](double t) {
        return A0 * exp((k1 - k2) * t);
    };

    std::vector<double> meanVal(nTimePoints);
    for (int i = 0; i < nTimePoints; i++) {
        meanVal[i] = compMean(times[i]);
    }

    title("k1 = 0.1, k2 = 0.01");
    xlabel("Time");
    ylabel("Amount of A");
    plot(times, a_amt, "-o");
    hold(on);
    plot(times, meanVal, "--r");

    ::matplot::legend({"Simulated A", "Calculated mean"});
    show();
    return 0;
}