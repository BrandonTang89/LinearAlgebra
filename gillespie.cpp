
#include <functional>
#include <random>
#include <vector>
#include <iostream>

#include "gillespie.h"
#include <cassert>

static constexpr double eps = 1e-12;
std::pair<std::vector<double>, std::vector<std::vector<double>>> glllespie(double t_final, std::vector<double> reactants, std::vector<Reaction> reactions, double volume) {
    std::random_device rd;  // obtain a random number from hardware
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    int nReactions = reactions.size();
    int nSpecies = reactants.size();
    std::vector<std::vector<double>> trajectory;
    std::vector<double> times;
    trajectory.push_back(reactants);
    times.push_back(0);

    double t = 0;
    while (t < t_final) {
        double r1 = dis(gen);
        double r2 = dis(gen);

        std::vector<double> cumulativePropensities(nReactions+1, 0);
        for (int i = 0; i < nReactions; i++) {
            cumulativePropensities[i+1] = reactions[i].propensity(reactants, volume);
            cumulativePropensities[i+1] += cumulativePropensities[i];
        }
        double a0 = cumulativePropensities[nReactions];
        if (a0 < eps) {
            // No more reactions possible
            break;
        }
        double tau = (1.0 / a0) * log(1.0 / r1);  // time to next reaction

        // we binary search for the index (idx) of the last reaction such that such that tau * a0 >= cumulativePropensities[idx]
        int idx = std::upper_bound(cumulativePropensities.begin(), cumulativePropensities.end(), r2 * a0) - 1 - cumulativePropensities.begin();
        assert(idx >= 0 && idx < nReactions);
        
        for (int i = 0; i < nSpecies; i++) {
            reactants[i] += reactions[idx].changes[i];
        }
        t += tau;
        trajectory.push_back(reactants);
        times.push_back(t);
    }

    return {times, trajectory};
}