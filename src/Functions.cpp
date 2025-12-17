#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>
#include <fstream>
#include <filesystem>
#include "Functions.hpp"
#include "Disk.h"
#include "Plan.h"
#include <iostream>
#include "Config.hpp"


uint32_t generate_seed() {
    std::random_device rd; // matériel / OS entropy si dispo
    uint32_t time_seed =
        static_cast<uint32_t>(std::chrono::high_resolution_clock::now()
            .time_since_epoch().count());

    // Combine rd et temps
    std::seed_seq seq { rd(), time_seed };
    std::vector<uint32_t> seeds(1);
    seq.generate(seeds.begin(), seeds.end());

    return seeds[0];
}

void place_grains(std::vector<Disk>& grains, const config::Config& cfg, const int task_index) {
    
    std::seed_seq seq{static_cast<uint32_t>(cfg.simulation.seed_value), static_cast<uint32_t>(task_index)};
    std::mt19937 gen(seq);
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    double x, y, vx, vy;
    int number_of_placed_grains = 0;
    int number_of_overlaps;

    // Global enforcement of p
    std::vector<bool> all_patches;
    if (cfg.patch.enforce_global_p) {
        // Total number of donor patches to assign
        const int total_patches = cfg.particles.N * cfg.patch.n;
        const int n_donor_patches = static_cast<int>(std::round(cfg.patch.p * total_patches));
        const int n_acceptor_patches = total_patches - n_donor_patches;

        // Create a vector with the appropriate number of donor (true) and acceptor (false) patches
        all_patches.insert(all_patches.end(), n_donor_patches, true);  // Donor patches
        all_patches.insert(all_patches.end(), n_acceptor_patches, false); // Acceptor patches

        // Shuffle the vector to randomize patch assignments across all grains
        std::shuffle(all_patches.begin(), all_patches.end(), gen);
    }


    // place grains
    while (number_of_placed_grains < cfg.particles.N)
    {
        number_of_overlaps = 0;
        double radius = cfg.particles.min_rad + dist(gen) * (cfg.particles.max_rad - cfg.particles.min_rad);
        
        x = -cfg.boundaries.lx / 2. + radius + dist(gen) * (cfg.boundaries.lx - 2. * radius);

        if (cfg.boundaries.oscillation.shape == 1 || cfg.boundaries.oscillation.shape == 2)
            y = -cfg.boundaries.ly/2. + cfg.boundaries.oscillation.correction + radius + dist(gen) * (cfg.boundaries.ly - 2. * radius - 2. * cfg.boundaries.oscillation.correction);
        else
            y = -cfg.boundaries.ly/2. + radius + dist(gen) * (cfg.boundaries.ly - 2. * radius);

        vx = 0.;
        vy = 0.;
        for (Disk &dsk : grains)
        {
            if (dsk.is_touching(x, y, radius))
            {
                number_of_overlaps++;
                break;
            }
        }
        if (number_of_overlaps == 0)
        {
            double mass = 4.0 / 3.0 * M_PI * radius * radius * radius * cfg.particles.density;
            std::vector<bool> patch_state(cfg.patch.n);
            if (cfg.patch.enforce_global_p) {
                
                // Assign patches to the current grain
                for (int k = 0; k < cfg.patch.n; k++) {
                    patch_state[k] = all_patches[number_of_placed_grains * cfg.patch.n + k];
                }
            } else {
                // Independent assignment per grain
                for (int k = 0; k < cfg.patch.n; k++) {
                    double rand_val = dist(gen);
                    if (rand_val < cfg.patch.p) {
                        patch_state[k] = true; // donneur
                    } else {
                        patch_state[k] = false; // accepteur
                    }
                }
            }
            double qmax = cfg.patch.q * cfg.patch.sat;
            grains.emplace_back(number_of_placed_grains, cfg.patch.n, cfg.patch.sat, qmax, patch_state, radius, mass, x, y, vx, vy);
            number_of_placed_grains++;
        }
    }

    if (cfg.output.save_layout)
    {
        std::ofstream myfile;
        std::string fileName = cfg.output.directory + "/initial_layouts/il_task_" + std::to_string(task_index) + ".txt";
        myfile.open(fileName);
        myfile.precision(10);
        for (Disk &dsk : grains)
        {
            myfile << dsk.index() << "\t"
                << dsk.r().x() << "\t" << dsk.r().y() << "\t"
                << dsk.v().x() << "\t" << dsk.v().y() << "\t"
                << dsk.theta() << "\t" << dsk.radius() << "\t";

            for (int k = 0; k < cfg.patch.n; k++) {
                myfile << dsk.getPatchStates().at(k) << "\t";
            }

            for (int k = 0; k < cfg.patch.n; k++) {
                myfile << dsk.getPatchCharges().at(k) << "\t";
            }

            myfile << std::endl;
        }
        myfile.close();
    }
}

// Pour ramener l'angle dans l'intervalle [0,2pi]
double norm_0_2pi(double angle) {
    angle = fmod(angle, 2.0*M_PI);
    return (angle < 0.0) ? angle + 2.0*M_PI : angle;
}

// Donne l'index du patch (0..n-1) d'un disque 
int patch_index_from_angle(double ang, int n) {
double angle = norm_0_2pi(ang);
double sector = 2.0*M_PI / (double)n;
int patch_index = (int)std::floor(angle / sector);
if (patch_index >= n) patch_index = n-1;
if (patch_index < 0)  patch_index = 0;
return patch_index;
}

double patch_angle_from_index(int patch_index, int n) {
    return ((patch_index+1) * (2.0*M_PI / (double)n) + patch_index * (2.0*M_PI / (double)n))/2.;
}


bool compute_disk_disk_contact(Disk& i_disk, Disk* j_disk_ptr, double i_deltan, Eigen::Vector3d& i_n, double i_kn, double i_e, double i_mu, bool toggle, double q, double qmax)
{
    bool ret = false;

    Disk* a = &i_disk;
    Disk* b = j_disk_ptr;

    //contact base
    i_n.normalize();
    Eigen::Vector3d v =  b->v() + b->radius()*b->w().cross(i_n) - (a->v() - a->radius()*a->w().cross(i_n));
    double vn_scalar = v.dot(i_n);
    Eigen::Vector3d vt = v - vn_scalar*i_n;
    double vt_scalar = vt.norm();
    Eigen::Vector3d t = (vt_scalar > 0.)? vt.normalized():Eigen::Vector3d(0.,1.,0.);
    
    //contact forces and torque
    double effectiveMass = (a->mass()*b->mass())/(a->mass()+b->mass());
    double eta = -2.*log(i_e)*sqrt(effectiveMass*i_kn/(log(i_e)*log(i_e)+M_PI*M_PI));
    
    //forces
    double Fn = i_kn*i_deltan + eta*vn_scalar;
    if(Fn > 0.)
    {
        a->add_force(Fn*i_n);
        b->add_force(-Fn*i_n);

        // =========================================================
        //    Détermination des PATCHES en contact (a et b)
        // =========================================================
        // Direction du point de contact côté a : -n
        // Direction du point de contact côté b : +n
        double ang_abs_b = std::atan2(i_n.y(), i_n.x());
        double ang_abs_a = norm_0_2pi(ang_abs_b + M_PI);

        // Passage dans le repère des particules
        double ang_a = norm_0_2pi(ang_abs_a - a->theta());
        double ang_b = norm_0_2pi(ang_abs_b - b->theta());

        int n_a = a->patch_count(); // (si n sont identiques, na==nb==n)
        int n_b = b->patch_count();
        
        int a_patch_index = patch_index_from_angle(ang_a, n_a);
        int b_patch_index = patch_index_from_angle(ang_b, n_b);

        // On a un transfert de charge si un patch donneur touche un patch accepteur (et réciproquement)

        bool a_patch = a->getPatchStates()[a_patch_index];
        bool b_patch = b->getPatchStates()[b_patch_index];



        //int dq_contact = -1; // !! Attention : mettre une charge postive casse les conditions if ci-dessous !! (car on vérifie pas la valeur absolue pour la saturation pour pas ralentir le programme)

        // DEBUG
        //std::cout << "Transfert de charge entre patch " << a_patch_index << " de la particule " << a->index() << " et patch " << b_patch_index << " de la particule " << b->index() << std::endl;
        //std::cout << "i_n: (" << i_n.x() << "," << i_n.y() << "Angles abs: a=" << ang_abs_a << " b=" << ang_abs_b << std::endl;
        

        
        
        a->addContact(a_patch_index, b->index(), b_patch_index, toggle);
        //b->addContact(b_patch_index, a->index(), a_patch_index, toggle); Je ne sais pas si c'est nécessaire de le faire des deux côtés


        if (a->getContactChargeTransfer(a_patch_index, b->index(), b_patch_index) == false && b->getContactChargeTransfer(b_patch_index, a->index(), a_patch_index) == false) // Si pas de transfert déjà effectué dans ce contact
        {
            //std::cout << "test" << std::endl;
            if (a_patch && !b_patch) // Si a est donneur et b accepteur
            {
                double a_charge = a->getPatchCharge(a_patch_index);
                double b_charge = b->getPatchCharge(b_patch_index);
                if (a_charge < qmax && -b_charge < qmax) // Si patch donneur/accepteur n'est pas saturé
                {
                    double charge = std::min(q, std::min(qmax - a_charge, qmax + b_charge));
                    a->add_charge(a_patch_index, charge);
                    b->add_charge(b_patch_index, -charge);
                    a->setContactChargeTransfer(a_patch_index, b->index(), b_patch_index, true);
                    b->setContactChargeTransfer(b_patch_index, a->index(), a_patch_index, true);
                    ret = true;

                }
            }
            else if ((!a_patch && b_patch)) // Si a est accepteur et b donneur
            {
                double a_charge = a->getPatchCharge(a_patch_index);
                double b_charge = b->getPatchCharge(b_patch_index);
                if (-a_charge < qmax && b_charge < qmax) // Si patch accepteur/donneur n'est pas saturé
                {
                    double charge = std::min(q, std::min(qmax + a_charge, qmax - b_charge));
                    a->add_charge(a_patch_index, -charge);
                    b->add_charge(b_patch_index, charge);
                    a->setContactChargeTransfer(a_patch_index, b->index(), b_patch_index, true);
                    b->setContactChargeTransfer(b_patch_index, a->index(), a_patch_index, true);
                    ret = true;
                }
            }
            // sinon : pas de transfert (couple inéligible)
        }
    }
    else
    {
        Fn=0.;
    }
    
    double Ft = 1000000.*vt_scalar;
    if(Ft > i_mu*Fn)
    {
        Ft = i_mu*Fn;
    }
    a->add_force(Ft*t);
    b->add_force(-Ft*t);
    
    //torque
    Eigen::Vector3d M = -Ft*a->radius()*i_n.cross(t);
    a->add_momentum(M);
    M = -Ft*b->radius()*i_n.cross(t);
    b->add_momentum(M);

    return ret;
}


void compute_disk_circular_piston_contact(Disk& i_disk, Disk* j_disk_ptr, double i_deltan, Eigen::Vector3d& i_n, double i_kn, double i_e, double i_mu, int grounded_walls)
{
    Disk* a = &i_disk;
    Disk* b = j_disk_ptr;
    
    //contact base
    i_n.normalize();
    Eigen::Vector3d v =  b->v() + b->radius()*b->w().cross(i_n) - (a->v() - a->radius()*a->w().cross(i_n));
    double vn_scalar = v.dot(i_n);
    Eigen::Vector3d vt = v - vn_scalar*i_n;
    double vt_scalar = vt.norm();
    Eigen::Vector3d t = (vt_scalar > 0.)? vt.normalized():Eigen::Vector3d(0.,1.,0.);
    
    //contact forces and torque
    double effectiveMass = (a->mass()*b->mass())/(a->mass()+b->mass());
    double eta = -2.*log(i_e)*sqrt(effectiveMass*i_kn/(log(i_e)*log(i_e)+M_PI*M_PI));
    
    //forces
    double Fn = i_kn*i_deltan + eta*vn_scalar;
    if(Fn > 0.)
    {
        a->add_force(Fn*i_n);
        b->add_force(-Fn*i_n);        
    }
    else
    {
        Fn=0.;
    }
    
    double Ft = 1000000.*vt_scalar;
    if(Ft > i_mu*Fn)
    {
        Ft = i_mu*Fn;
    }
    a->add_force(Ft*t);
    b->add_force(-Ft*t);
    
    //torque
    Eigen::Vector3d M = -Ft*a->radius()*i_n.cross(t);
    a->add_momentum(M);
    M = -Ft*b->radius()*i_n.cross(t);
    b->add_momentum(M);

    
    if (grounded_walls == 1)
    {
        double angle = std::atan2(i_n.y(), i_n.x()) + M_PI; // Direction du point de contact côté disque
        angle = norm_0_2pi(angle - i_disk.theta()); // Passage dans le repère de la particule
        int n = i_disk.patch_count();
        int patch_index = patch_index_from_angle(angle, n);
        i_disk.reset_patch_charge(patch_index);
    }
    else if (grounded_walls == 2)
    {
        i_disk.reset_patch_charges();
    }
}


void compute_disk_wall_contact(Disk& i_disk, double i_deltan, Eigen::Vector3d& i_n, double i_kn, double i_e, double i_mu, int grounded_walls)
{
    //contact base
    i_n.normalize();
    Eigen::Vector3d v = -i_disk.v() + i_disk.radius()*i_disk.w().cross(i_n);
    double vn_scalar = v.dot(i_n);
    Eigen::Vector3d vt = v - vn_scalar*i_n;
    double vt_scalar = vt.norm();
    Eigen::Vector3d t = (vt_scalar > 0.)? vt.normalized():Eigen::Vector3d(0.,1.,0.);
    
    //contact forces and torque
    double effectiveMass = i_disk.mass();
    double eta = -2.*log(i_e)*sqrt(effectiveMass*i_kn/(log(i_e)*log(i_e)+M_PI*M_PI));
    
    //forces
    double Fn = i_kn*i_deltan + eta*vn_scalar;
    if(Fn > 0.)
    {
        i_disk.add_force(Fn*i_n);
    }
    else
    {
        Fn=0.;
    }
    
    double Ft = 1000000.*vt_scalar;
    if(Ft > i_mu*Fn)
    {
        Ft = i_mu*Fn;
    }
    i_disk.add_force(Ft*t);
    
    //torque
    Eigen::Vector3d M = -Ft*i_disk.radius()*i_n.cross(t);
    i_disk.add_momentum(M);

    if (grounded_walls == 1)
    {
        double angle = std::atan2(i_n.y(), i_n.x()) + M_PI; // Direction du point de contact côté disque
        angle = norm_0_2pi(angle - i_disk.theta()); // Passage dans le repère de la particule
        int n = i_disk.patch_count();
        int patch_index = patch_index_from_angle(angle, n);
        i_disk.reset_patch_charge(patch_index);
    }
    else if (grounded_walls == 2)
    {
        i_disk.reset_patch_charges();
    }
}

void compute_screened_coulomb_interaction(Disk& a, Disk& b, double ke, double inv_lambda, double r_soft, double r_cut2, std::vector<double>& patch_angle_cache, int n_patch, double fast_r_cut2)
{
    const double rad_a = a.radius();
    const double rad_b = b.radius();
    const Eigen::Vector3d n = b.r() - a.r();
    if (n.squaredNorm() > fast_r_cut2) return;

    const double theta_a = a.theta();
    const double theta_b = b.theta();

    const Eigen::Vector3d ar = a.r();
    const Eigen::Vector3d br = b.r();
    Eigen::Vector3d F_tot = Eigen::Vector3d::Zero();
    Eigen::Vector3d M_a_tot = Eigen::Vector3d::Zero();
    Eigen::Vector3d M_b_tot = Eigen::Vector3d::Zero();

    std::vector<Eigen::Vector3d> r_a(n_patch), r_b(n_patch), pos_a(n_patch), pos_b(n_patch);
    std::vector<double> q_a(n_patch), q_b(n_patch);


    // precompute patch positions & charges
    for (int i = 0; i < n_patch; ++i) {
        double sin_a, cos_a, sin_b, cos_b;
        sincos(theta_a + patch_angle_cache[i], &sin_a, &cos_a);
        sincos(theta_b + patch_angle_cache[i], &sin_b, &cos_b);
        r_a[i] = Eigen::Vector3d{rad_a * cos_a, rad_a * sin_a, 0.0};
        r_b[i] = Eigen::Vector3d{rad_b * cos_b, rad_b * sin_b, 0.0};
        pos_a[i] = ar + r_a[i];
        pos_b[i] = br + r_b[i];
        q_a[i] = a.getPatchCharge(i);
        q_b[i] = b.getPatchCharge(i);
  
    }

    for (int i = 0; i < n_patch; i++) {
        const double qi = q_a[i];
        if (qi == 0.0) continue;
        const Eigen::Vector3d& posAi = pos_a[i];
        const Eigen::Vector3d& rAi   = r_a[i];
        const double ke_qi = ke * qi;

        for (int j = 0; j < n_patch; j++) {
            const double qj = q_b[j];
            if (qj == 0.0) continue;

            Eigen::Vector3d rvec = posAi - pos_b[j];
            const double rvec2 = rvec.squaredNorm();
            if (rvec2 > r_cut2) continue;
            double r = std::sqrt(rvec2);

            // adoucissement pour éviter r → 0 (le contact gère déjà le court-portée)
            double r_eff = std::max(r, r_soft);
            double inv_reff = 1.0 / r_eff;

            double expf = std::exp(-r_eff * inv_lambda);

            // Force de Yukawa: F = k q_i q_j e^{-r/λ} ( 1/r^2 + 1/(λ r) ) r̂
            double mag =  ke_qi * qj * expf * (inv_reff * inv_reff + inv_lambda * inv_reff);

            Eigen::Vector3d F = (mag/r) * rvec;
            F_tot += F;

            // Moments (couples) autour des centres
            // tau = r x F (en 2D: composante z)
            M_a_tot += rAi.cross(F);      // sur A
            M_b_tot -= r_b[j].cross(F);            
        }
    }
    a.add_force( F_tot);
    b.add_force(-F_tot);
    a.add_momentum(M_a_tot);
    b.add_momentum(M_b_tot);
}

double Ep_compute_screened_coulomb_interaction(Disk& a, Disk& b, double ke, double inv_lambda, double r_soft, double r_cut2, std::vector<double>& patch_angle_cache, int n_patch, double fast_r_cut2)
{
    double Ep = 0.0;
    const double rad_a = a.radius();
    const double rad_b = b.radius();
    const Eigen::Vector3d n = b.r() - a.r();
    if (n.squaredNorm() > fast_r_cut2) return 0.0;

    const double theta_a = a.theta();
    const double theta_b = b.theta();

    const Eigen::Vector3d ar = a.r();
    const Eigen::Vector3d br = b.r();
    Eigen::Vector3d F_tot = Eigen::Vector3d::Zero();
    Eigen::Vector3d M_a_tot = Eigen::Vector3d::Zero();
    Eigen::Vector3d M_b_tot = Eigen::Vector3d::Zero();

    std::vector<Eigen::Vector3d> r_a(n_patch), r_b(n_patch), pos_a(n_patch), pos_b(n_patch);
    std::vector<double> q_a(n_patch), q_b(n_patch);


    // precompute patch positions & charges
    for (int i = 0; i < n_patch; ++i) {
        double sin_a, cos_a, sin_b, cos_b;
        sincos(theta_a + patch_angle_cache[i], &sin_a, &cos_a);
        sincos(theta_b + patch_angle_cache[i], &sin_b, &cos_b);
        r_a[i] = Eigen::Vector3d{rad_a * cos_a, rad_a * sin_a, 0.0};
        r_b[i] = Eigen::Vector3d{rad_b * cos_b, rad_b * sin_b, 0.0};
        pos_a[i] = ar + r_a[i];
        pos_b[i] = br + r_b[i];
        q_a[i] = a.getPatchCharge(i);
        q_b[i] = b.getPatchCharge(i);
  
    }

    for (int i = 0; i < n_patch; i++) {
        const double qi = q_a[i];
        if (qi == 0.0) continue;
        const Eigen::Vector3d& posAi = pos_a[i];
        const Eigen::Vector3d& rAi   = r_a[i];
        const double ke_qi = ke * qi;

        for (int j = 0; j < n_patch; j++) {
            const double qj = q_b[j];
            if (qj == 0.0) continue;

            Eigen::Vector3d rvec = posAi - pos_b[j];
            const double rvec2 = rvec.squaredNorm();
            if (rvec2 > r_cut2) continue;
            double r = std::sqrt(rvec2);

            // adoucissement pour éviter r → 0 (le contact gère déjà le court-portée)
            double r_eff = std::max(r, r_soft);
            double inv_reff = 1.0 / r_eff;

            double k_expf = ke_qi * qj * std::exp(-r_eff * inv_lambda);

            // Force de Yukawa: F = k q_i q_j e^{-r/λ} ( 1/r^2 + 1/(λ r) ) r̂
            double mag =  k_expf * (inv_reff * inv_reff + inv_lambda * inv_reff);

            Eigen::Vector3d F = (mag/r) * rvec;
            F_tot += F;

            Ep += k_expf * inv_reff;

            // Moments (couples) autour des centres
            // tau = r x F (en 2D: composante z)
            M_a_tot += rAi.cross(F);      // sur A
            M_b_tot -= r_b[j].cross(F);
        }
    }
    a.add_force( F_tot);
    b.add_force(-F_tot);
    a.add_momentum(M_a_tot);
    b.add_momentum(M_b_tot);

    return Ep;
}


double kinetic_energy(Disk& dsk)
{
    double KE_trans = 0.5 * dsk.mass() * dsk.v().squaredNorm();
    double KE_rot = 0.5 * dsk.inertia() * dsk.w().squaredNorm();
    return KE_trans + KE_rot;
}