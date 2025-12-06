# DEM_ChargeBuildUp

A Discrete Element Method (DEM) simulation tool for charge buildup analysis.

## Overview

This project focuses on Discrete Element Method (DEM) simulations for charge buildup analysis.

## Requirements

This project requires the following libraries:

- **Eigen**: Linear algebra library for C++
- **toml++**: TOML configuration file parser for C++
- **Voro++**: Voronoi tessellation library
  - ⚠️ **Important**: Voro++ must be compiled before use. Follow the compilation instructions in the Voro++ documentation.

## Installation

```bash
# Install Eigen (example for Ubuntu/Debian)
sudo apt-get install libeigen3-dev

# Install toml++
# Follow instructions from: https://github.com/marzer/tomlplusplus

# Compile Voro++
# Download from: https://math.lbl.gov/voro++/
cd voro++-0.4.6
make
sudo make install
```

## Usage

```bash
# Example with 4 processes
mpiexec -n 4 ./dem.exe
```



## Contact

**Louis Cleen**  
Email: louis.cleen@gmail.com
