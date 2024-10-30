#ifndef HELPERS_H
#define HELPERS_H

#include <cmath>

constexpr double EPS = 1e-9;
constexpr double INF = 1e9;
inline bool isCloseTo(double a, double b) {
    return std::fabs(a - b) < EPS;
}

#endif  // HELPERS_H