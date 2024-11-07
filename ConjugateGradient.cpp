#include <functional>
#include <iostream>

#include "ConjugateGradient.h"
#include "Matrix.h"
#include "Vec.h"

Vec conjugateGradient(const Matrix& A, const Vec& b, int maxIter, double tol) {
    int n = b.N;
    Vec x(n, 0.0);  // Initial guess is a zero vector
    Vec r = b;      // r = b - Ax (initially 0)
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

Vec conjugateGradientPreconditioned(const Matrix& A, const Vec& b, int maxIter, double tol, std::function<Vec(Vec)> multMInv) {
    // Precondition: multMInv(v) = M_inv * v where M is a positive definite preconditioning matrix.
    int n = b.N;
    Vec x(n, 0.0);  // Initial guess is a zero vector
    Vec rh = b;
    Vec ph = rh;
    double rhMinvrh = rh.dot(multMInv(rh)); // r_T M_-1 r

    for (int i = 0; i < maxIter && rh.norm() > tol; ++i) {
        Vec Aph = A * ph;
        double alpha = rhMinvrh / (ph.dot(Aph));
        x = x + alpha * ph;
        rh = rh - alpha * Aph;
        double rhMinvrh_new = rh.dot(multMInv(rh));
        double beta = rhMinvrh_new / rhMinvrh;
        ph = multMInv(rh) + beta * ph;
        rhMinvrh = rhMinvrh_new;
    }

    return x;
}

Vec conjugateGradientPreconditionedDiagonal(const Matrix& A, const Vec& b, int maxIter, double tol){
    std::vector<double> diagA = A.diag();
    auto multMInv = [&diagA](Vec v) {
        Vec result(v.N);
        for (int i = 0; i < v.N; i++) {
            result[i] = v[i] / diagA[i]; // diagA is non-zero since A is positive definite
        }
        return result;
    };

    Vec v(b.N, 1.0);
    return conjugateGradientPreconditioned(A, b, maxIter, tol, multMInv);
}