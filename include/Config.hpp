#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <toml.hpp>
#include <Eigen/Dense>
#include <string>
#include <array>
#include <vector>

namespace config {


struct SimulationConfig
{
    std::string name;
    std::string seed;
    std::string seed_source;
    uint32_t seed_value;
    int n_runs;
    int current_run;
    double threshold_cage;
};

struct TimeConfig
{
    double total;
    double rec_start;
    double mean_time;
    int fps;
    double dt;
    int n_snapshots;
};

struct ParticlesConfig
{
    std::vector<int> N_list;
    int N;
    std::string initial_layout;
    double min_rad;
    double max_rad;
    double density;
    std::vector<double> g;
};

struct PatchConfig
{
    double p;
    bool enforce_global_p;
    std::vector<int> sat_list;
    int sat;
    int n;
    std::vector<double> q_list;
    double q;
    int grounded_walls;
};

struct BoundariesConfig
{
    double lx;
    double ly;

    struct OscillationConfig
    {
        int  shape; // 0 = straight, 1 = inclined, 2 = circular
        double correction;
        double angle;
        double A_bot;
        double A_top;
        double f_bot;
        double f_top;
    } oscillation;
};

struct ContactConfig
{
    double e;
    double mu;
    double kn;
};

struct ElectrostaticConfig
{
    double eps_r;
    double lambda_e;
};

struct OutputConfig
{
    std::string directory;
    std::string save_mode_str;
    int save_mode;
    bool save_layout;
};

struct Config
{
    SimulationConfig    simulation;
    TimeConfig          time;
    ParticlesConfig     particles;
    PatchConfig         patch;
    BoundariesConfig    boundaries;
    ContactConfig       contact;
    ElectrostaticConfig electrostatic;
    OutputConfig        output;
};



Config from_toml(const toml::table& tbl, std::string seed = "", int rank = 0); // Lit config.toml et remplit une structure Config.
Config load_config(const std::string& content, std::string seed = "", int rank = 0);
Config load_config_file(const std::string& path, std::string seed = "", int rank = 0);
std::vector<double> get_double_values(const toml::node_view<const toml::node>& node, int rank = 0);
std::vector<int> get_int_values(const toml::node_view<const toml::node>& node, int rank = 0);
void get_parameters(Config& cfg, int task_index); // Récupère les paramètres pour un run


} // namespace config

#endif
