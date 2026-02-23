#define _USE_MATH_DEFINES
#include "Config.hpp"
#include "Functions.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <toml.hpp>
#include <vector>

namespace config {

Config from_toml(const toml::table &tbl, std::string seed, std::string output_path, int rank) {
    Config cfg;

    auto sim = tbl["simulation"];
    cfg.simulation.name = sim["name"].value_or("Unknown");
    cfg.simulation.seed = sim["seed"].value_or("");

    if (!seed.empty()) {
        cfg.simulation.seed_value = static_cast<uint32_t>(std::stoul(seed));
        cfg.simulation.seed_source = "from command line";
    } else if (!cfg.simulation.seed.empty()) {
        cfg.simulation.seed_value = static_cast<uint32_t>(std::stoul(cfg.simulation.seed));
        cfg.simulation.seed_source = "from config file";
    } else {
        if (rank == 0)
            cfg.simulation.seed_value = generate_seed();
        cfg.simulation.seed_source = "generated";
    }

    cfg.simulation.n_runs = sim["n_runs"].value_or(1);
    cfg.simulation.threshold_cage = sim["threshold_cage"].value_or(0.4);

    auto time = tbl["time"];
    cfg.time.total = time["total"].value_or(1.0);
    cfg.time.rec_start = time["rec_start"].value_or(0.0);
    cfg.time.mean_time = time["mean_time"].value_or(-1.0);
    if (cfg.time.rec_start + cfg.time.mean_time > cfg.time.total) {
        throw std::runtime_error(
            "mean_time exceeds the available simulation time after rec_start.");
    } else if (cfg.time.mean_time < 0.0) {
        cfg.time.mean_time = cfg.time.total - cfg.time.rec_start;
    }
    cfg.time.fps = time["fps"].value_or(100);
    cfg.time.dt = time["dt"].value_or(5e-05);
    cfg.time.n_snapshots = std::floor((cfg.time.total - cfg.time.rec_start) * cfg.time.fps);

    auto particles = tbl["particles"];
    cfg.particles.N_list = get_int_values(particles["N"], rank);
    cfg.particles.detect_transitions = particles["detect_transitions"].value_or(false);
    cfg.particles.number_of_values_around_transition =
        particles["number_of_values_around_transition"].value_or(10);
    cfg.particles.initial_layout = particles["initial_layout"].value_or("random");
    // cfg.particles.min_rad = particles["min_rad"].value_or(0.0005);
    // cfg.particles.max_rad = particles["max_rad"].value_or(0.0005);
    cfg.particles.radii = get_double_values(particles["radii"], rank);
    cfg.particles.density = particles["density"].value_or(10000.0);
    auto g_array = particles["g"].as_array();
    if (g_array != nullptr && g_array->size() == 3) {
        cfg.particles.g.clear();
        for (const auto &val : *g_array) {
            cfg.particles.g.push_back(val.value_or(0.0));
        }
    } else {
        cfg.particles.g = {0.0, 0.0, 0.0};
    }

    auto patch = tbl["patch"];
    cfg.patch.p = patch["p"].value_or(0.5);
    cfg.patch.enforce_global_p = patch["enforce_global_p"].value_or(true);
    cfg.patch.sat_list = get_int_values(patch["sat"], rank);
    cfg.patch.n = patch["n"].value_or(7);
    cfg.patch.q_list = get_double_values(patch["q"], rank);
    cfg.patch.grounded_walls = patch["grounded_walls"].value_or(0);
    cfg.patch.tau = patch["tau"].value_or(0.0);
    cfg.patch.tau_enabled = (cfg.patch.tau > 0.0);
    cfg.patch.tau_leak = patch["tau_leak"].value_or(0.0);
    cfg.patch.tau_leak_enabled = (cfg.patch.tau_leak > 0.0);

    if (cfg.patch.tau_enabled && cfg.time.dt > cfg.patch.tau / 10.0) {
        if (rank == 0) {
            std::cerr << "Warning: Time step dt is large compared to tau. Consider reducing dt for "
                         "better stability."
                      << std::endl;
        }
    }
    if (cfg.patch.tau_leak_enabled && cfg.time.dt > cfg.patch.tau_leak / 10.0) {
        if (rank == 0) {
            std::cerr << "Warning: Time step dt is large compared to tau_leak. Consider reducing "
                         "dt for better stability."
                      << std::endl;
        }
    }

    auto boundaries = tbl["boundaries"];
    cfg.boundaries.lx = boundaries["lx"].value_or(0.03);
    cfg.boundaries.ly = boundaries["ly"].value_or(0.06);
    cfg.boundaries.use_normalized_length = boundaries["use_normalized_length"].value_or(false);
    cfg.boundaries.normalized_length = get_double_values(boundaries["normalized_length"], rank);

    auto oscillation = boundaries["oscillation"];
    cfg.boundaries.oscillation.shape = oscillation["shape"].value_or(0);
    cfg.boundaries.oscillation.angle = oscillation["angle"].value_or(0.0);
    cfg.boundaries.oscillation.A_bot = oscillation["A_bot"].value_or(0.005);
    cfg.boundaries.oscillation.A_top = oscillation["A_top"].value_or(0.005);
    cfg.boundaries.oscillation.f_bot = oscillation["f_bot"].value_or(10.0);
    cfg.boundaries.oscillation.f_top = oscillation["f_top"].value_or(10.0);

    cfg.boundaries.oscillation.correction = 0.0;
    if (cfg.boundaries.oscillation.shape == 1)
        cfg.boundaries.oscillation.correction =
            std::tan(cfg.boundaries.oscillation.angle) *
            (cfg.boundaries.lx / 2.); // hauteur du triangle formé par le plan incliné
    else if (cfg.boundaries.oscillation.shape == 2)
        cfg.boundaries.oscillation.correction =
            cfg.boundaries.lx -
            (cfg.boundaries.lx / 2.) * std::sqrt(3.0); // correction de la hauteur de la boite

    auto contact = tbl["contact"];
    cfg.contact.e = contact["e"].value_or(0.9);
    cfg.contact.mu = contact["mu"].value_or(0.06);
    cfg.contact.kn = contact["kn"].value_or(10000);
    cfg.contact.delta_max_over_R = contact["delta_max_over_R"].value_or(0.01);

    auto electrostatic = tbl["electrostatic"];
    cfg.electrostatic.eps_r = electrostatic["eps_r"].value_or(1.0);
    cfg.electrostatic.lambda_e = electrostatic["lambda_e"].value_or(5e-4);

    auto output = tbl["output"];
    cfg.output.directory = output_path + output["directory"].value_or("data");
    cfg.output.save_mode_str = output["save_mode"].value_or("reduced");
    if (cfg.output.save_mode_str == "full")
        cfg.output.save_mode = 2;
    else if (cfg.output.save_mode_str == "reduced")
        cfg.output.save_mode = 1;
    else
        cfg.output.save_mode = 0;

    cfg.output.save_initial_layout = output["save_initial_layout"].value_or(false);
    cfg.output.save_final_layout = output["save_final_layout"].value_or(false);
    return cfg;
}

Config load_config(const std::string &content, std::string seed, std::string output_path,
                   int rank) {
    toml::table tbl;
    try {
        tbl = toml::parse(content);
    } catch (const toml::parse_error &err) {
        if (rank == 0) {
            std::cerr << "Error parsing TOML content: " << err.description() << "\n"
                      << "At " << err.source().begin << std::endl;
        }
        throw;
    }

    return from_toml(tbl, seed, output_path, rank);
}

Config load_config_file(const std::string &file_path, std::string seed, std::string output_path,
                        int rank) {
    toml::table tbl;
    try {
        tbl = toml::parse_file(file_path);
    } catch (const toml::parse_error &err) {
        if (rank == 0) {
            std::cerr << "Error parsing TOML file: " << err.description() << "\n"
                      << "At " << err.source().begin << std::endl;
        }
        throw;
    }

    return from_toml(tbl, seed, output_path, rank);
}

std::vector<double> get_double_values(const toml::node_view<const toml::node> &node, int rank) {
    std::vector<double> values;

    // cas 1 : scalaire
    if (node.is_floating_point() || node.is_integer()) {
        values.push_back(node.value<double>().value());
    }

    // cas 2 : liste
    else if (node.is_array()) {
        for (auto &v : *node.as_array()) {
            values.push_back(v.value<double>().value());
        }
    }

    // cas 3 : table {start, end, n}
    else if (node.is_table()) {
        auto t = *node.as_table();

        double start = t["start"].value<double>().value();
        double end = t["end"].value<double>().value();
        int n = t["n"].value<int>().value();

        if (n <= 1) {
            if (rank == 0) {
                std::cerr << "Warning: n should be greater than 1 for range specification. "
                             "Returning start value only."
                          << std::endl;
            }
            return std::vector<double>{start};
        }

        values.reserve(n);
        for (int i = 0; i < n; i++) {
            values.push_back(start + (end - start) * i / (n - 1));
        }
    }

    return values;
}

std::vector<int> get_int_values(const toml::node_view<const toml::node> &node, int rank) {
    std::vector<int> values;

    // cas 1 : scalaire
    if (node.is_integer()) {
        values.push_back(node.value<int>().value());
    }

    // cas 2 : liste
    else if (node.is_array()) {
        for (auto &v : *node.as_array()) {
            values.push_back(v.value<int>().value());
        }
    }

    // cas 3 : table {start, end, n}
    else if (node.is_table()) {
        auto t = *node.as_table();

        int start = t["start"].value<int>().value();
        int end = t["end"].value<int>().value();
        int n = t["n"].value<int>().value();

        if (n <= 1) {
            if (rank == 0) {
                std::cerr << "Warning: n should be greater than 1 for range specification. "
                             "Returning start value only."
                          << std::endl;
            }
            return std::vector<int>{start};
        }

        values.reserve(n);
        for (int i = 0; i < n; i++) {
            if ((end - start) * i % (n - 1) != 0)
                if (rank == 0)
                    std::cerr << "Warning: non-integer step size in int range specification. "
                                 "Values will be truncated."
                              << std::endl;

            values.push_back(start + (end - start) * i / (n - 1));
        }
    }

    return values;
}

int determine_transition_N(const Config &cfg, int N_index) {
    double sweep_coeff =
        0.3; // Coefficient pour déterminer la plage autour de la transition à simuler
    double a = -1275.17229151;
    double b = 200.8994774;

    double A_top = cfg.boundaries.oscillation.A_top;
    double A_bot = cfg.boundaries.oscillation.A_bot;
    double ly = cfg.boundaries.ly;
    double lx = cfg.boundaries.lx;

    if (cfg.particles.min_rad != cfg.particles.max_rad) {
        std::cerr
            << "Warning: determine_transition_N is currently implemented for monodisperse systems. "
               "min_rad and max_rad are different. Using min_rad for transition determination."
            << std::endl;
    }
    double r = cfg.particles.min_rad;

    double S_tot = lx * (ly + A_bot + A_top);
    double Sg = M_PI * r * r;

    double normalized_length = ly / r;
    double pf = (normalized_length - b) / a; // Estimation du pf à partir de la longueur normalisée
    double N_transition = (pf * S_tot) / Sg; // Estimation de N à la transition

    int N_min = static_cast<int>(std::floor(N_transition * (1.0 - sweep_coeff)));
    int N_max = static_cast<int>(std::ceil(N_transition * (1.0 + sweep_coeff)));

    if (N_min < 20) {
        return -1; // Indique que la transition n'est pas détectable dans une plage raisonnable
    }

    std::vector<int> N_values =
        linspace_int(N_min, N_max, cfg.particles.number_of_values_around_transition);

    return N_values[N_index];
}

void get_parameters(Config &cfg, int task_index) {
    if (task_index < 0)
        throw std::invalid_argument("Task index must be non-negative.");

    // Get parameters for patch.sat and patch.q
    int n_R = cfg.particles.radii.size();
    int n_L;
    if (cfg.boundaries.use_normalized_length)
        n_L = cfg.boundaries.normalized_length.size();
    else
        n_L = 1; // If not using normalized length, we consider it as a single fixed parameter

    int n_N;
    if (cfg.particles.detect_transitions) {
        // If detect_transitions is enabled, we ignore the specified N_list and determine N based on
        // the transition values
        n_N = cfg.particles.number_of_values_around_transition; // We will simulate a number of
                                                                // values around the transition N
    } else {
        n_N = cfg.particles.N_list.size();
    }
    int n_sat = cfg.patch.sat_list.size();
    int n_q = cfg.patch.q_list.size();
    int n_run = cfg.simulation.n_runs;

    int total_combinations = n_R * n_L * n_N * n_sat * n_q * n_run;

    if (n_R == 0 || n_L == 0 || n_sat == 0 || n_q == 0 || n_N == 0)
        throw std::runtime_error("At least one parameter list for particles.radii, "
                                 "boundaries.normalized_length, particles.N, "
                                 "patch.sat, or patch.q is empty.");

    if (task_index >= total_combinations)
        throw std::out_of_range("Task index exceeds total number of parameter combinations.");

    int R_index = task_index % n_R;
    task_index /= n_R;
    int L_index = task_index % n_L;
    task_index /= n_L;
    int N_index = task_index % n_N;
    task_index /= n_N;
    int sat_index = task_index % n_sat;
    task_index /= n_sat;
    int q_index = task_index % n_q;
    task_index /= n_q;
    int run_index = task_index % n_run;
    task_index /= n_run;

    // Assign selected parameters to cfg
    cfg.particles.min_rad = cfg.particles.radii[R_index];
    cfg.particles.max_rad = cfg.particles.radii[R_index];
    if (cfg.boundaries.use_normalized_length)
        cfg.boundaries.ly = cfg.boundaries.normalized_length[L_index] * (cfg.particles.max_rad);
    if (cfg.particles.detect_transitions) {
        // If detect_transitions is enabled, we ignore the specified N_list and determine N based on
        // the transition values
        int transition_N = determine_transition_N(
            cfg, N_index); // This function should be implemented to return the
        //  transition N based on preliminary simulations
        cfg.particles.N = transition_N; // You can also consider adding some variability around
                                        // the transition N if desired
    } else {
        cfg.particles.N = cfg.particles.N_list[N_index];
    }

    cfg.patch.sat = cfg.patch.sat_list[sat_index];
    cfg.patch.q = cfg.patch.q_list[q_index];
    cfg.simulation.current_run = run_index;

    compute_kn(cfg);
    compute_dt(cfg);

    if (cfg.output.save_mode_str == "hybrid") {
        if (run_index == 0)
            cfg.output.save_mode = 2; // full
        else
            cfg.output.save_mode = 1; // reduced
    }
}

} // namespace config