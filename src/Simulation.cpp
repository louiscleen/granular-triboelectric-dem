#define _USE_MATH_DEFINES
#include "Simulation.hpp"
#include "Config.hpp"
#include "Functions.hpp"
#include "Disk.h"
#include "Plan.h"
#include "Cell.h"
#include "SimulationIO.hpp"

#include <mpi.h>
#include <voro++.hh>

#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <random>
#include <filesystem>
#include <fstream>


AggregatedResult simulation(const config::Config& cfg, const int task_index) {
    // To avoid dereferencing, define frequently used parameters:

    // Electrostatic constants
    const double eps0 = 8.854187817e-12;
    const double lambda_e = cfg.electrostatic.lambda_e;
    const double k_e = 1.0 / (4.0 * M_PI * eps0 * cfg.electrostatic.eps_r);
    const double r_cut_e  = 4.0 * lambda_e;   // cutoff pour efficacité
    const double r_soft_e = 0.25 * ((cfg.particles.min_rad + cfg.particles.max_rad) * 0.5); // adoucissement numérique
    const double q = cfg.patch.q;

    const double threshold_cage = cfg.simulation.threshold_cage;
    
    // Contact parameters
    const double kn = cfg.contact.kn;
    const double e = cfg.contact.e;
    const double mu = cfg.contact.mu;

    // Time parameters
    const double total_time = cfg.time.total;
    const double fps = cfg.time.fps;
    const double start_time = cfg.time.rec_start;
    const double dt = cfg.time.dt;


    // Container dimensions
    double lx = cfg.boundaries.lx;
    double ly = cfg.boundaries.ly;
    bool circular = false;

    // Oscillation parameters
    const double A_bot = cfg.boundaries.oscillation.A_bot;
    const double y0_bot = -ly / 2.;                                          // position moyenne piston inférieur
    const double omega_bot = 2. * M_PI * cfg.boundaries.oscillation.f_bot;   // fréquence angulaire piston inférieur
    const double A_top = cfg.boundaries.oscillation.A_top;
    const double y0_top = ly / 2.;                                           // position moyenne piston supérieur
    const double omega_top = 2. * M_PI * cfg.boundaries.oscillation.f_top;   // fréquence angulaire piston supérieur 

    // Gravity
    Eigen::Vector3d gravity = Eigen::Vector3d(cfg.particles.g[0], cfg.particles.g[1], cfg.particles.g[2]);

    // Output parameters
    int save_mode = cfg.output.save_mode;

    std::unique_ptr<SimulationWriter> writer; // Simulation writer pointer (for reduced data or full data)

    // More variables
    int cellIndex;
    Eigen::Vector3d r, n;
    bool toggle = true; // Optimization toggle for contact calculations
    int frame_id = 0;   // Frame index for output

    // Results storage
    std::vector<Result> results;
    results.reserve(cfg.time.n_snapshots);


    // Boundary walls initialization
    std::vector<Plan> walls;

    if (cfg.boundaries.oscillation.shape == 1) {
        walls.emplace_back(0., -ly / 2., -1*std::sin(cfg.boundaries.oscillation.angle), std::cos(cfg.boundaries.oscillation.angle)); // bottom wall
        walls.emplace_back(0., ly / 2., std::sin(cfg.boundaries.oscillation.angle), -1.*std::cos(cfg.boundaries.oscillation.angle));  // top wall
        walls.emplace_back(-lx / 2., 0., 1., 0.);   // left wall
        walls.emplace_back(lx / 2., 0., -1., 0.);   // right wall

    } else if (cfg.boundaries.oscillation.shape == 2) {
        walls.emplace_back(-lx / 2., 0., 1., 0.);   // left wall
        walls.emplace_back(lx / 2., 0., -1., 0.);   // right wall
        // top and bottom walls are replaced by disks in circular case

    } else {
        walls.emplace_back(0., -ly / 2., 0, 1);      // bottom wall
        walls.emplace_back(0., ly / 2., 0., -1.);    // top wall
        walls.emplace_back(-lx / 2., 0., 1., 0.);    // left wall
        walls.emplace_back(lx / 2., 0., -1., 0.);    // right wall
    }

    // Disks for circular oscillation
    double position_bottom_disk = 0.0;
    double position_top_disk = 0.0;
    Disk bottom_disk;
    Disk top_disk;
    if (cfg.boundaries.oscillation.shape ==  2) {
        circular = true;
        position_bottom_disk = -ly/2. - (lx/2.)*std::sqrt(3.0);
        position_top_disk = ly/2. + (lx/2.)*std::sqrt(3.0);
        bottom_disk = Disk(-2, 0, 0, std::vector<bool>(), lx, 1.0, 0.0, position_bottom_disk, 0.0, 0.0);
        top_disk = Disk(-1, 0, 0, std::vector<bool>(), lx, 1.0, 0.0, position_top_disk, 0.0, 0.0);
    }



    // Initial placement of grains
    std::vector<Disk> grains;
    grains.reserve(cfg.particles.N);
    place_grains(grains, cfg, task_index); // task_index to ensure different placements across tasks (for RNG seeding)

    
    // Correct container height for oscillating walls
    if (cfg.boundaries.oscillation.shape == 1) {
        ly += A_bot + cfg.boundaries.oscillation.A_top + 2*cfg.boundaries.oscillation.correction; 
    } else if (cfg.boundaries.oscillation.shape == 2) {
        ly += A_bot + 2*lx + A_top;
    } else {
        ly += A_bot + A_top; 
    }

    // Initialize Voro++ container
    double nx = cbrt(cfg.particles.N) * (lx / ly); // Nombre de cellules en x pour le conteneur Voro++
    double ny = cbrt(cfg.particles.N) * (ly / lx); // Nombre de cellules en y pour le conteneur Voro++
    int init_mem_block = 15; // Mémoire initiale pour Voro++
    voro::container_poly con(-lx/2., lx/2., -ly/2., ly/2., 0.0, 1.0, nx, ny, 1,false,false,false, init_mem_block); // conteneur Voro++


    // Linked cells initialization
    double cellSize = 2.2 * cfg.particles.max_rad;
    int extendedNeighbor_radius = 3;

    int nCellx = lx / cellSize;
    int nCelly = ly / cellSize;
    int nCell = nCellx * nCelly;
    std::vector<Cell> cellules;
    cellules.reserve(nCell);
    double dx = lx / nCellx;
    double dy = ly / nCelly;

    int ix, iy, jx, jy;
    for (int i = 0; i < nCell; i++)
    {
        iy = i / nCellx;
        ix = i % nCellx;
        cellules.emplace_back(i);
        for (int j = 0; j < i; j++)
        {
            jy = j / nCellx;
            jx = j % nCellx;
            int delta_x = abs(ix - jx);
            int delta_y = abs(iy - jy);

            if (delta_x < 2 && delta_y < 2)
            {
                cellules[i].add_neighbor(&cellules[j]);
            }
            else if (delta_x < extendedNeighbor_radius && delta_y < extendedNeighbor_radius)
            {
                cellules[i].add_neighbor_of_neighbor(&cellules[j]);
            }
        }
    }

    

    // Reduced & full data output initialization
    std::string destination_full_dir;
    if (save_mode == 2 ) {
        destination_full_dir = cfg.output.directory + "/full_data/task_" + std::to_string(task_index);
        std::filesystem::create_directories(destination_full_dir);

    } else if (save_mode == 1) {
        writer = std::make_unique<SimulationWriter>(cfg.output.directory + "/reduced_data/out_task_" + std::to_string(task_index) + ".bin", cfg);
    }

    // Minimal data output initialization
    std::ofstream file_minimal_data;
    file_minimal_data.open(cfg.output.directory + "/minimal_data/md_task_" + std::to_string(task_index) + ".txt");
    file_minimal_data.precision(10); 
    file_minimal_data << "Time" << "\t" << "N_c" << "\t" << "KE" << "\t" << "charge_transfers" << std::endl;

    int charge_transfers = 0;


    // Gravity force initialization
    for (Disk &dsk : grains) {
        dsk.add_gravity_force(gravity);
    }

    // ================================================================
    // ===================== Main simulation loop =====================
    // ================================================================
    for (double time = 0.; time < total_time; time += dt)
    {
        if (circular) {
            // ------------ Circular plan ------------
            bottom_disk.set_position(0.0, y0_bot - (lx/2.)*std::sqrt(3.0) + A_bot * std::sin(omega_bot * time));
            top_disk.set_position(0.0, y0_top + (lx/2.)*std::sqrt(3.0) + A_top * std::sin(omega_top * time + M_PI));
        } else {
            // ------------ Straight and inclined plan ------------
            walls[0].set_position(0.0, y0_bot + A_bot * std::sin(omega_bot * time)); // bottom wall
            walls[1].set_position(0.0, y0_top + A_top * std::sin(omega_top * time + M_PI)); // top wall
        }


        // linked cells
        for (Cell &cl : cellules)
        {
            cl.set_head_of_list(nullptr);
        }

        // grains
        for (Disk &dsk : grains)
        {
            // update velocities
            dsk.update_velocity(0.5 * dt);

            // update positions
            dsk.update_position(dt);

            // class in linked cells
            cellIndex = (int)((dsk.r().x() + lx / 2.) / dx) + (int)((dsk.r().y() + ly / 2.) / dy) * nCellx;
            Cell *cl_ptr = &cellules[cellIndex];

            dsk.set_linked_cell(cl_ptr);
            dsk.set_linked_disk(cl_ptr->head_of_list());
            cl_ptr->set_head_of_list(&dsk);

            // reset force
            dsk.reset_force();

            // set gravity force
            dsk.add_gravity_force(gravity);
        }

        // contact detection
        for (Disk &dsk : grains)
        {
            // in same cell
            Disk *other_dsk_ptr = dsk.linked_cell()->head_of_list();
            while (other_dsk_ptr != nullptr)
            {
                if (other_dsk_ptr->index() < dsk.index())
                {
                    n = dsk.r() - other_dsk_ptr->r();
                    double delta = (dsk.radius() + other_dsk_ptr->radius()) - n.norm();

                    add_screened_coulomb(dsk, *other_dsk_ptr, k_e, lambda_e, r_soft_e, r_cut_e, q);
                    
                    if (delta > 0.) {
                        charge_transfers += compute_contact(dsk, other_dsk_ptr, delta, n, kn, e, mu, toggle);
                    }
                }
                other_dsk_ptr = other_dsk_ptr->linked_disk();
            }

            // in neighbor cells
            for (Cell *cl : dsk.linked_cell()->neighbors())
            {
                Disk *other_dsk_ptr = cl->head_of_list();
                while (other_dsk_ptr != nullptr)
                {
                    n = dsk.r() - other_dsk_ptr->r();
                    double delta = (dsk.radius() + other_dsk_ptr->radius()) - n.norm();

                    add_screened_coulomb(dsk, *other_dsk_ptr, k_e, lambda_e, r_soft_e, r_cut_e, q);

                    if (delta > 0.)
                    {
                        charge_transfers += compute_contact(dsk, other_dsk_ptr, delta, n, kn, e, mu, toggle);
                    }
                    other_dsk_ptr = other_dsk_ptr->linked_disk();
                }
            }

            // in neighbor of neighbor cells
            for (Cell *cl : dsk.linked_cell()->neighbors_of_neighbors())
            {
                Disk *other_dsk_ptr = cl->head_of_list();
                while (other_dsk_ptr != nullptr)
                {
                    add_screened_coulomb(dsk, *other_dsk_ptr, k_e, lambda_e, r_soft_e, r_cut_e, q);
                    other_dsk_ptr = other_dsk_ptr->linked_disk();
                }
            }

            // with walls
            for (Plan &pl : walls)
            {
                n = pl.n();
                r = dsk.r() - pl.r();
                double delta = dsk.radius() - r.dot(n);

                if (delta > 0.)
                {
                    compute_contact(dsk, delta, n, kn, e, mu);
                }
            }

            if (circular)
            {
                n = dsk.r() - bottom_disk.r();
                double delta = (dsk.radius() + bottom_disk.radius()) - n.norm();

                if (delta > 0.) {
                    compute_contact(dsk, &bottom_disk, delta, n, kn, e, mu);
                }

                n = dsk.r() - top_disk.r();
                delta = (dsk.radius() + top_disk.radius()) - n.norm();

                if (delta > 0.) {
                    compute_contact(dsk, &top_disk, delta, n, kn, e, mu);
                }
            }
        }

        // update velocity and position
        for (Disk &dsk : grains)
        {
            dsk.update_velocity(0.5 * dt);

            // clean contacts
            for (auto it = dsk.getContacts().begin(); it != dsk.getContacts().end();)
            {
                if (it->state != toggle)
                {
                    it = dsk.getContacts().erase(it);
                }
                else
                {
                    ++it;
                }
            }
        }


        // Recording
        if (time >= frame_id / fps) {
            if (time >= start_time) {
                
                int partic_valide = 0;
                int partic_caged = 0;

                double KE = 0.;
                
                for (Disk &dsk : grains)
                {
                    dsk.setCaged(false); // reset caged state
                    con.put(dsk.index(), dsk.r().x(), dsk.r().y(), 0.5, dsk.radius());  // insert particle into Voro++ container
                    KE += kinetic_energy(dsk);
                }
                voro::c_loop_all cl(con);   
                voro::voronoicell cell;    

                
                if (cl.start()) do {
                    if (con.compute_cell(cell, cl)) {
                        int id = cl.pid();
                        //double vg = 4.0/3.0 * M_PI * std::pow(grains[id].radius(), 3);
                        partic_valide++;
                        double sg = M_PI * grains[id].radius()*grains[id].radius();
                        if (sg/cell.volume() > threshold_cage) {
                            grains[id].setCaged(true);
                            partic_caged++;
                        }
                    }
                } while (cl.inc());

                file_minimal_data << time << "\t" << partic_caged << "\t" << KE <<  "\t" << charge_transfers << "\n";
                
                results.emplace_back(partic_caged, KE); // Using C++20 aggregate initialization

                con.clear(); // clear Voro++ container for next use

                
                if (save_mode == 2)
                    writer = std::make_unique<SimulationWriter>(destination_full_dir + "/out_" + std::to_string(frame_id) + ".txt", cfg);

                if (save_mode != 0) {
                    // Prepare particle data for reduced output
                    std::vector<ParticleData> particles;
                    particles.reserve(grains.size());
                    for (Disk& dsk : grains) {
                        particles.emplace_back(dsk.index(), dsk.r().x(), dsk.r().y(), dsk.v().x(), dsk.v().y(), dsk.theta(), dsk.radius(), dsk.isCaged()); // Using C++20 aggregate initialization
                    }
                    writer->write_snapshot(particles, time);
                }

            }
            frame_id++;
        } // end record

        toggle = !toggle;
    } // end main simulation loop
    file_minimal_data.close();



    // Compute mean values over specified time interval
    if (cfg.time.mean_time == 0.0)
        return AggregatedResult{static_cast<double>(results.back().caged_particles), results.back().kinetic_energy};

    double partic_caged_mean = 0;
    double KE_mean = 0.0;
    int mean_frames = static_cast<int>(cfg.time.mean_time * fps); // Number of frame corresponding to the time interval

    for (std::size_t i = results.size() - mean_frames; i < results.size(); ++i)
    {
        partic_caged_mean += results[i].caged_particles;
        KE_mean += results[i].kinetic_energy;
    }
    return AggregatedResult{partic_caged_mean / mean_frames, KE_mean / mean_frames};
}










