#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <Eigen/Dense>
#include <vector>

#include "Config.hpp"

class Cell;
class Disk;
class Plan;

unsigned int generate_seed();
void place_grains(std::vector<Disk> &grains, const config::Config &cfg,
                  const int task_index); // task_index to ensure different placements across tasks
                                         // (for RNG seeding)
void save_layout(std::vector<Disk> &grains, const config::Config &cfg, const int task_index,
                 const std::string &prefix);

bool compute_disk_disk_contact(Disk &, Disk *, double, Eigen::Vector3d &, double, double, double,
                               bool, double, double);
void compute_disk_circular_piston_contact(Disk &, Disk *, double, Eigen::Vector3d &, double, double,
                                          double, int);
void compute_disk_wall_contact(Disk &, double, Eigen::Vector3d &, double, double, double, int);

void compute_screened_coulomb_interaction(Disk &, Disk &, double, double, double, double,
                                          std::vector<double> &, int, double);
double Ep_compute_screened_coulomb_interaction(Disk &, Disk &, double, double, double, double,
                                               std::vector<double> &, int, double);
double kinetic_energy(Disk &);

std::vector<int> linspace_int(int a, int b, int n);

void compute_kn(config::Config &cfg);
void compute_dt(config::Config &cfg);

#endif // FUNCTIONS_HPP