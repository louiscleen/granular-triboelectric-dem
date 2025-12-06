#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "Config.hpp"

struct Result {
    int caged_particles;
    double kinetic_energy;
};

struct AggregatedResult {
    double mean_caged_particles;
    double mean_kinetic_energy;
};

AggregatedResult simulation(const config::Config& cfg, const int task_index);


#endif // SIMULATION_HPP
