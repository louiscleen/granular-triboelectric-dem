#ifndef SIMULATIONIO_HPP
#define SIMULATIONIO_HPP

#include "Config.hpp"
#include "Constants.hpp"

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>
// #include <stdexcept>

// ===========================================================================
// STRUCTURES DE DONNÉES
// ===========================================================================

// Struct pour les particules
struct ParticleData {
    int index;
    double x, y;
    double vx, vy;
    double theta;
    double radius;
    double charge_global;
    double charge_abs;
    bool isCaged;
};

// Paramètres de la simulation
struct SimulationParams {
    uint32_t seed;
    double total_time;
    double start_time;
    double fps;

    int N;
    double min_rad;
    double max_rad;
    double density;
    double p;
    int sat;
    int n;
    double q;
    char shape; // 's' = straight, 'i' = inclined, 'c' = circular
    double angle;
    double A_bot;
    double A_top;
    double f_bot;
    double f_top;
};

// Header global du fichier
struct FileHeader {
    char magic[8] = {'D', 'C', 'B', 'U'};
    uint32_t format_version = 1; // version du format de fichier
    char simulation_name[65];    // nom de la simulation
    char version[9];             // version du logiciel
    uint32_t n_snapshots;        // nombre total de snapshots
    uint32_t param_size;         // taille de param
};

// Header d'un snapshot
struct SnapshotHeader {
    double time; // temps du snapshot
};

// ===========================================================================
// OBJET : SimulationWriter  → écriture au fil de l’eau
// ===========================================================================
class SimulationWriter {
  public:
    SimulationWriter(const std::string &filename, const config::Config &cfg)
        : save_mode_(cfg.output.save_mode) {
        if (save_mode_ == 2) {
            out.open(filename);
        } else if (save_mode_ == 1) {
            out.open(filename, std::ios::binary);
            if (!out)
                throw std::runtime_error("Impossible d'ouvrir le fichier");

            // Préparer le header global
            FileHeader h;
            h.n_snapshots = cfg.time.n_snapshots;
            h.param_size = sizeof(SimulationParams);

            // Copie de la version du logiciel dans un tableau fixe
            std::memset(h.version, 0, 9);
            std::strncpy(h.version, std::string(VERSION).c_str(), 8);

            // Copie du nom de simulation dans un tableau fixe
            std::memset(h.simulation_name, 0, 65);
            std::strncpy(h.simulation_name, cfg.simulation.name.c_str(), 64);

            // Écriture du header
            out.write((char *)&h, sizeof(h));

            // Préparer les paramètres de la simulation
            SimulationParams params;
            params.seed = cfg.simulation.seed_value;
            params.total_time = cfg.time.total;
            params.start_time = cfg.time.rec_start;
            params.fps = cfg.time.fps;
            params.N = cfg.particles.N;
            params.min_rad = cfg.particles.min_rad;
            params.max_rad = cfg.particles.max_rad;
            params.density = cfg.particles.density;
            params.p = cfg.patch.p;
            params.sat = cfg.patch.sat;
            params.n = cfg.patch.n;
            params.q = cfg.patch.q;
            params.shape =
                cfg.boundaries.oscillation.shape; // 's' = straight, 'i' = inclined, 'c' = circular
            params.angle = cfg.boundaries.oscillation.angle;
            params.A_bot = cfg.boundaries.oscillation.A_bot;
            params.A_top = cfg.boundaries.oscillation.A_top;
            params.f_bot = cfg.boundaries.oscillation.f_bot;
            params.f_top = cfg.boundaries.oscillation.f_top;

            // Écriture des paramètres
            out.write((char *)&params, sizeof(params));
        }
    }

    // Écrit un snapshot (header + particules)
    void write_snapshot(const std::vector<ParticleData> &particles, double time) {
        if (save_mode_ == 2) {
            // Écriture en mode texte complet (non binaire)
            out.precision(10);
            for (const auto &p : particles) {
                out << p.index << "\t" << p.x << "\t" << p.y << "\t" << p.vx << "\t" << p.vy << "\t"
                    << p.theta << "\t" << p.radius << "\t" << p.charge_global << "\t"
                    << p.charge_abs << "\t" << p.isCaged << "\n";
            }
            return;
        } else if (save_mode_ == 1) {
            // Écriture en mode binaire réduit
            SnapshotHeader sh; // Header du snapshot

            sh.time = time;

            out.write((char *)&sh, sizeof(sh));

            // Données des particules
            out.write((char *)particles.data(), particles.size() * sizeof(ParticleData));
        }
    }

    ~SimulationWriter() {
        if (out.is_open())
            out.close();
    }

  private:
    std::ofstream out;
    int save_mode_ = 1;
};

// ===========================================================================
// OBJET : SimulationReader → lecture séquentielle du fichier
// ===========================================================================
class SimulationReader {
  public:
    SimulationReader(const std::string &filename, int save_mode = 1) : save_mode_(save_mode) {

        if (save_mode_ == 2) {
            in.open(filename);
            if (!in)
                throw std::runtime_error("Impossible d'ouvrir le fichier");

        } else if (save_mode_ == 1) {
            in.open(filename, std::ios::binary);
            if (!in)
                throw std::runtime_error("Impossible d'ouvrir le fichier");

            // Lire header
            in.read((char *)&header, sizeof(header));
            if (std::string(header.magic, 4) != "DCBU")
                throw std::runtime_error("Format de fichier invalide");
            if (header.format_version != 1)
                throw std::runtime_error("Version de format non supportée");

            // Lire paramètres
            if (header.param_size != sizeof(SimulationParams))
                throw std::runtime_error("Taille des paramètres inattendue");
            in.read((char *)&params, sizeof(params));

            N_particles_ = params.N;
        }
    }

    const FileHeader &getHeader() const { return header; }
    const SimulationParams &getParams() const { return params; }

    bool read_snapshot(std::vector<ParticleData> &particles, double &time) {
        if (save_mode_ == 2) {
            // Lecture en mode texte complet (non binaire)
            particles.clear();
            std::string line;
            while (std::getline(in, line)) {
                ParticleData p;
                std::istringstream iss(line);
                if (!(iss >> p.index >> p.x >> p.y >> p.vx >> p.vy >> p.theta >> p.radius >>
                      p.charge_global >> p.charge_abs >> p.isCaged)) {
                    break; // fin de fichier ou erreur
                }
                particles.push_back(p);
            }
            time = -1.; // not stored in text mode... could be improved

            return !particles.empty();

        } else if (save_mode_ == 1) {
            // Lecture en mode binaire
            if (in.peek() == EOF)
                return false;

            SnapshotHeader sh;
            in.read((char *)&sh, sizeof(sh));
            if (!in)
                return false;

            particles.resize(N_particles_);
            in.read((char *)particles.data(), N_particles_ * sizeof(ParticleData));
            if (!in)
                return false;

            time = sh.time;
            return true;
        }
        return false;
    }

    ~SimulationReader() {
        if (in.is_open())
            in.close();
    }

  private:
    std::ifstream in;
    FileHeader header;
    SimulationParams params;
    int save_mode_ = 1;
    int N_particles_;
};

#endif // SIMULATIONIO_HPP
