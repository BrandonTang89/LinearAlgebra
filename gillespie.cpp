
#include <functional>
#include <random>
#include <vector>

#include "gillespie.h"

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

        double a0 = 0;
        std::vector<double> cumulativePropensities(nReactions, 0);
        for (int i = 0; i < nReactions; i++) {
            cumulativePropensities[i] += reactions[i].propensity(reactants, volume);
            a0 += cumulativePropensities[i];
            if (i > 0) {
                cumulativePropensities[i] += cumulativePropensities[i - 1];
            }
        }
        double tau = (1.0 / a0) * log(1.0 / r1);  // time to next reaction

        // we binary search for the index (idx) of the first reaction such that such that tau * a0 >= cumulativePropensities[idx]
        int idx = std::lower_bound(cumulativePropensities.begin(), cumulativePropensities.end(), r2 * a0) - cumulativePropensities.begin();
        for (int i = 0; i < nSpecies; i++) {
            reactants[i] += reactions[idx].changes[i];
        }
        t += tau;
        trajectory.push_back(reactants);
        times.push_back(t);
    }

    return {times, trajectory};
}