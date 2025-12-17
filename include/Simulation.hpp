#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "Config.hpp"

struct Result {
    int caged_particles;
    double kinetic_energy;
    double electrostatic_energy;
};

struct AggregatedResult {
    double mean_caged_particles;
    double mean_kinetic_energy;
    double mean_electrostatic_energy;
    double charge_transfers_normalized;
    double charge_abs_total;

};

AggregatedResult simulation(const config::Config& cfg, const int task_index);


#endif // SIMULATION_HPP
