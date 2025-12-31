#include "Config.hpp"
#include "Constants.hpp"
#include "DefaultConfig.hpp"
#include <cxxopts.hpp>
#include <filesystem>
#include <iostream>
#include <mpi.h>
#include <string>

void parsing_options(int argc, char *argv[], std::string &config_file, std::string &seed,
                     std::string &output_path, int &exit_requested, int rank) {
    cxxopts::Options options("DEM Charge build-up",
                             "A Discrete Element Method simulation for charge build-up.");

    options.add_options()("h,help", "Display this information.")("v, version", "Display version")(
        "c,config", "Config file (default: config.toml)",
        cxxopts::value<std::string>()->default_value(std::string(DEFAULT_CONFIG_FILE)))(
        "s,seed", "Global random seed for all MPI processes",
        cxxopts::value<unsigned int>())("o, output",
                                        "Output path (default: current working directory) where "
                                        "'directory' from config file will be created",
                                        cxxopts::value<std::string>()->default_value(""))
        //("e, ETA", "Estimated time of arrival based on number of MPI processes")
        ;

    auto parsedArgs = options.parse(argc, argv);

    config_file = parsedArgs["config"].as<std::string>();

    if (parsedArgs.count("seed")) {
        seed = std::to_string(parsedArgs["seed"].as<unsigned int>());
    } else {
        seed = "";
    }

    if (rank == 0) {
        if (parsedArgs.count("help")) {
            std::cout << options.help() << std::endl;
            exit_requested = 1;
            return;
        }

        if (parsedArgs.count("version")) {
            std::cout << "DEM_ChargeBuildUp version " << VERSION << std::endl;
            exit_requested = 1;
            return;
        }

        if (config_file == DEFAULT_CONFIG_FILE &&
            !std::filesystem::exists(std::filesystem::path(DEFAULT_CONFIG_FILE))) {
            std::cout << "No config.toml file found. Generating default file..." << std::endl;

            config::write_default_config();

            std::cout << "You can now modify the file and restart the program." << std::endl;
            exit_requested = 1;
            return;
        } else if (!std::filesystem::exists(config_file)) {
            std::cerr << "Config file not found: " << config_file << std::endl;
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
    }

    output_path = parsedArgs["output"].as<std::string>();
    if (!output_path.empty() && output_path.back() != '/' && output_path.back() != '\"' &&
        output_path.back() != '\\') {
        output_path += '/';
    }

    // if (parsedArgs.count("ETA")) {
    //     std::cout << "Estimated time of arrival: " << std::endl;
    //     MPI_Finalize();
    // }
}
