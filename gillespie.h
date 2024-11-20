#pragma once

#include <functional>
#include <random>
#include <vector>

struct Reaction {
    std::function<double(std::vector<double>, double)> propensity;
    // probability of reaction happening in [t, t+dt_) = propensity(reactants, volume) * dt
    std::vector<int> changes;  // changes in reactants
};

std::pair<std::vector<double>, std::vector<std::vector<double>>> glllespie(double t_final, std::vector<double> reactants, std::vector<Reaction> reactions, double volume = 1.0);