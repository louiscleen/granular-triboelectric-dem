#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <Eigen/Dense>
#include <vector>

#include "Config.hpp"


class Cell;
class Disk;
class Plan;

unsigned int generate_seed();
void place_grains(std::vector<Disk>& grains, const config::Config& cfg, const int task_index); // task_index to ensure different placements across tasks (for RNG seeding)
void compute_contact(Disk&,Disk*,double,Eigen::Vector3d&,double,double,double, bool);
void compute_contact(Disk&,Disk*,double,Eigen::Vector3d&,double,double,double);
void compute_contact(Disk&,double,Eigen::Vector3d&,double,double,double);

void add_screened_coulomb(Disk&, Disk&, double, double, double, double, double);
double kinetic_energy(Disk&);


#endif // FUNCTIONS_HPP