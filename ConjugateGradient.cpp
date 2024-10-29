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