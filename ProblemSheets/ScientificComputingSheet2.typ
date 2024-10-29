#import "@preview/codelst:2.0.1": sourcecode
#let title = "Scientific Computing Sheet 2"
#let author ="Brandon Tang"
#set document(title: title, author: author)
#set par(justify: true)
#set page(numbering: "1/1", number-align: right,)
#let tut(x) = [#block(x, stroke: blue, radius: 1em, inset: 1.5em, width: 100%)]
#let pblock(x) = [#block(x, stroke: rgb("#e6c5fc") + 0.03em, fill: rgb("#fbf5ff"), radius: 0.3em, inset: 1.5em, width: 100%)]
#let gblock(x) = [#block(x, stroke: rgb("#5eb575") + 0.03em, fill: rgb("#e3fae9"), radius: 0.3em, inset: 1.5em, width: 100%)]
#align(center)[
  #block(text(weight: 700, 1.75em, title))
  #v(1em, weak: true)
  #text(weight: 550, 1.1em, author)
]

#sourcecode([```cpp
#include "ConjugateGradient.h"
#include "Vec.h"
#include "Matrix.h"
#include <iostream>

Vec conjugateGradient(const Matrix& A, const Vec& b, int maxIter, double tol){
    int n = b.N;
    Vec x(n, 0.0); // Initial guess is a zero vector
    Vec r = b - A * x;
    Vec p = r;

    for (int i = 0; i < maxIter && r.norm() > tol; ++i) {
        Vec Ap = A * p;
        double alpha = r.dot(r) / p.dot(Ap);
        x = x + alpha * p;
        Vec r_new = r - alpha * Ap;
        double beta = r_new.dot(r_new) / r.dot(r);
        p = r_new + beta * p;
        r = r_new;
    }

    return x;
}
```])

#sourcecode([```cpp
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

    return 0;
}
```])

#pblock([
Solution to Ax = b using Conjugate Gradient: -0.147829 0.155398 0.0584987 -0.025608 0.00839462 0.0177763 0.0018769 -0.00137001 0.00396464 0.00210073 

Residual: 3.82899e-16
])