#ifndef DEFAULTCONFIG_HPP
#define DEFAULTCONFIG_HPP

#include "Constants.hpp"
#include <string>
#include <fstream>
#include <iostream>
#include <filesystem>

namespace config
{
inline std::string generate_default_config_toml()
{
    return R"(#########################################
#         Simulation DEM Config
#########################################
[simulation]
name = "DEM simulation" # 64 characters max
seed = ""               # Seed (empty for random generation), overriden by --seed command line argument if provided. Always saved in output files.
n_runs = 1 				# Number of times a simulation is repeated
param_sweep = "grid"	# How to explore the parameter space ("grid" for cartesian product) # To be added: LHS, corner case testing
threshold_cage = 0.4    # Threshold value of the packing fraction for considering a particle as caged in 2D (0.4 means 40% of the Voronoi cell is occupied by the particle)


#########################################
#             Parameters
#########################################
[time]
total = 10.0    		# Total simulation time 		[s]
rec_start = 0.0			# Capture start time			[s]
mean_time = -1.0        # Time over which the mean values are computed [s] 
# Note that mean values are computed over the interval [total - mean_time, total]
# If set to 0.0, only the final values at time = total are considered, if set to < 0.0, mean values are computed over [rec_start to total]
fps = 100 				# Framerate 					[frames/s]
dt = 1e-05 				# Time step						[s]
# Note that the number of snapshots is computed as: n_snapshots = (total-rec_start) * fps

[particles]
N = 50 					# Number of particles
# N can be set as a list or table to run multiple simulation with the following syntax:
# N = [50, 100, 120] is valid
# N = { start = 1, end = 100, n = 100 } is also valid

# Initial positions: "random" or "load" (to load from initial_layout.txt)
initial_layout = "random"	
min_rad = 0.0005 		# Minimum radius 				[m]
max_rad = 0.0005 		# Maximum radius 				[m]
density = 2500.0 		# Density				    	[kg/m³]
# For bronze, use density = 9500 kg/m³, for glass density = 2500 kg/m³, for m≈0.0005kg density = 954929.0 kg/m³
g = [0.0, 0.0, 0.0]		# Gravity						[m/s²]

[patch]
p = 0.5 				# Probability of donor patches
enforce_global_p = 1    # 0 means independent patch assignment per grain, 1 means global p is enforced (total number of donor patches across all grains is approx p*N*n)
# Note that if set to 0, the actual fraction of donor patches may differ from p due to statistical fluctuations (especially for small N)
sat = 10				# Maximum charge per patch
# sat can be set as a list or table to run multiple simulation with the following syntax:
# sat = [5, 10, 12] is valid
# sat = { start = 0, end = 10, n = 6 } is also valid

n = 7 					# Number of patches per grain
q = 2e-11 		    	# Elementary charge per patch 	[C]
# q can be set as a list or table to run multiple simulation with the following syntax:
# q = [2e-11, 4e-11, 8e-11] is valid
# q = { start = 0.1, end = 1.0, n = 10 } is also valid
tau = 0.0            # Global relaxation time on the particle for charge redistribution among patches (set to 0 to disable) [s]
tau_leak = 0.0          # Relaxation time for patch charge leakage to the outside (set to 0.0 to disable leakage) [s]
grounded_walls = 0      # If set to 0, walls are insulating (no charge transfer between particles and walls)
                        # If set to 1, walls are considered grounded (particles patches can transfer charge to walls, which remain at zero potential)
                        # If set to 2, walls are considered grounded and particles lose all their patch charges upon contact with a wall
                        

[boundaries] 			# 
lx = 0.03 				# Container width				[m]
ly = 0.06 				# Container	height				[m]
[boundaries.oscillation]
shape = 0          		# Piston shape, 0 = straight, 1 = inclined, 2 = circular)
angle = 0.0 			# Inclination angle (only for inclined shape)	[rad]
A_bot = 0.005 			# Bottom wall amplitude 		[m]  
A_top = 0.005 			# Top wall amplitude 			[m]		
f_bot = 10 				# Bottom wall frequency 		[Hz] 			
f_top = 10 				# Top wall frequency 			[Hz]			

[contact]
e = 0.9 				# Restitution coefficient 
mu = 0.6				# Friction coefficient 
kn = 1000.0 			# Normal stiffness  			[N/m]

[electrostatic]
eps_r = 1.0             # Relative permittivity
lambda_e = 5e-4 		# Screening length      		[m]
# Note that r_cut_e (cut-off radius) is defined as 4.0 * lambda_e                   [m]
# Note that r_soft_e (digital softening) is defined as 0.25 * (min_rad+max_rad)/2   [m]



#########################################
#             Output
#########################################
[output]
directory = "data"      # Data directory name
                        # The output directory will be created inside the current working directory. To specify a different path, use --output command line argument.
                        # Example: --output $GLOBALSCRATCH/$SLURM_JOB_ID
                        # Note that if the folder already exists, the program will stop to avoid overwriting existing data (the directory is always created by the program)
                        
# Controls how much simulation data is written to disk.
# - "full":    save particle positions, velocities, states and auxiliary data in seperate .txt files for each frame (Useful for debugging by reproducing a simulation from a seed)
# - "hybrid":  if n_runs=1, same as "full", if n_runs > 1, save full data for the first run and reduced data for the others
# - "reduced": same as above but in an single binary file per simulation (all snapshots in one file)
# - "minimal": save only the analysis in one global file & one minimal file per simulation (no particle data) 
# Note that full & reduced data also include the minimal data.
save_mode = "full"      
save_layout = false  	# Save initial layout (with information on particle patches) to a file "initial_layout_<task_index>.txt"

# Binary file format (used in "reduced" save_mode):
# For each simulation, all snapshots are stored in a single binary file with the following structure:

# +------------------------------+
# | Global file header           |  
# +------------------------------+
# | Snapshot 0 header            |
# | Snapshot 0 data              |
# +------------------------------+
# | Snapshot 1 header            |
# | Snapshot 1 data              |
# +------------------------------+
# | ...                          |)";
}

inline void write_default_config()
{
    std::string default_config(DEFAULT_CONFIG_FILE);
    std::ofstream file(default_config);
    if (!file)
    {
        std::cerr << "Unable to create the file: " << default_config << std::endl;
        return;
    }

    file << generate_default_config_toml();
    std::cout << "Configuration file generated: " << default_config << std::endl;
}
} // namespace config

#endif