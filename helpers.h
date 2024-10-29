#ifndef HELPERS_H
#define HELPERS_H

#include <cmath>

const double EPS = 1e-9;
inline bool isCloseTo(double a, double b) {
    return std::fabs(a - b) < EPS;
}

#endif // HELPERS_H