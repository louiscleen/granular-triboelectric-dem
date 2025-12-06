#include <toml.hpp>
#include <mpi.h>

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <filesystem>
#include <chrono>

#include "Constants.hpp"
#include "CommandLineParser.hpp"
#include "Config.hpp"
#include "Simulation.hpp"
#include "Functions.hpp"


struct TaskInfo
{
    int index;
    int rank;
    double duration;
};

struct RankInfo
{
    int rank;
    int num_tasks;
    double total_duration;
    double idle_time;
};





int main(int argc, char* argv[])
{
    // Initializing MPI
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);       // rank --> id du processus [0, size-1]
    MPI_Comm_size(MPI_COMM_WORLD, &size);       // size --> nombre de processus
    int exit_requested = 0;

    if (rank == 0)
        std::cout << "=== DEM Charge build-up with MPI ===" << std::endl;

    std::string config_file = "";
    std::string temp_seed = "";  // Should not be used, only for CommandLineParser 
    parsing_options(argc, argv, config_file, temp_seed, exit_requested, rank);

    // Broadcast exit_requested to all processes
    MPI_Bcast(&exit_requested, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (exit_requested) {
        MPI_Finalize();
        return 0;
    }


    std::string config_content;

    if (rank == 0) {
        std::ifstream f(config_file);
        config_content.assign((std::istreambuf_iterator<char>(f)),
                    std::istreambuf_iterator<char>());
    }

    int config_content_size = config_content.size();
    MPI_Bcast(&config_content_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    config_content.resize(config_content_size);

    MPI_Bcast(config_content.data(), config_content_size, MPI_CHAR, 0, MPI_COMM_WORLD);

    config::Config cfg;

    try {
        cfg = config::load_config(config_content, temp_seed, rank);
    }
    catch (const toml::parse_error& err) {
        if (rank == 0) {
            std::cerr << "Config loading failed: " << config_file << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // Broadcast the global seed to all processes
    MPI_Bcast(&cfg.simulation.seed_value, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);





    // =========================================================================
    //  Master-Worker Task Allocation
    //  (Simple version, master does not perform simulations)
    // =========================================================================

    if (rank == 0) {

        if (size < 2) {
            std::cerr << "This program requires at least 2 MPI processes (one master and at least one worker)." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        } 

        auto start_time = std::chrono::high_resolution_clock::now();
        
        std::cout << "Job name: " << cfg.simulation.name << std::endl;
        std::cout << "# Number of processes: " << size << " (1 master + " << size - 1 << " workers)." << std::endl;
        std::cout << "# Config loaded from: " << config_file << std::endl;
        std::cout << "# Using global seed: " << cfg.simulation.seed_value << " (" << cfg.simulation.seed_source << ")" << std::endl;

        // ----- MASTER -----
        int num_tasks = cfg.particles.N_list.size() * cfg.patch.sat_list.size() * cfg.patch.q_list.size() * cfg.simulation.n_runs;
        int num_tasks_one_run = cfg.particles.N_list.size() * cfg.patch.sat_list.size() * cfg.patch.q_list.size();
        std::cout << "# Total number of tasks: " << num_tasks << " (" << num_tasks_one_run << " per run)." << "\n" << std::endl;

        // Copier le fichier de configuration dans le dossier de sortie
        if (std::filesystem::exists(cfg.output.directory)) {
            std::cerr << "Error: the folder: \"" << cfg.output.directory << "\" already exists." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        try {
            std::filesystem::create_directories(cfg.output.directory); // si le dossier n'existe pas
            std::filesystem::copy_file(config_file, cfg.output.directory + "/" + config_file);
        }
        catch (const std::filesystem::filesystem_error& e) {
            std::cerr << "Error: " << e.what() << std::endl;
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        std::string destination_minimal_dir = cfg.output.directory + "/minimal_data";
        std::filesystem::create_directories(destination_minimal_dir);

        if (cfg.output.save_mode == 1) {
            std::string destination_reduced_dir = cfg.output.directory + "/reduced_data";
            std::filesystem::create_directories(destination_reduced_dir);
        } else if (cfg.output.save_mode == 2) {
            std::string destination_full_dir = cfg.output.directory + "/full_data";
            std::filesystem::create_directories(destination_full_dir);
        }

        if (cfg.output.save_layout) {
            std::string destination_layout_dir = cfg.output.directory + "/initial_layouts";
            std::filesystem::create_directories(destination_layout_dir);
        }


        int next_task = 0;
        std::vector<int> task_indices(num_tasks, -1);
        std::vector<AggregatedResult> results(num_tasks, {-1., -1.}); // Pour contenir moyennes des résultats
        std::vector<TaskInfo> durations(num_tasks, {-1, -1, 0.0});

        MPI_Status status;

        // envoyer une tâche initiale à chaque worker
        for (int worker = 1; worker < size && next_task < num_tasks; ++worker) {
            MPI_Send(&next_task, 1, MPI_INT, worker, TAG_WORK, MPI_COMM_WORLD); // envoyer l'indice de la tâche
            next_task++;
        }
        


        int tasks_completed = 0;

        // boucle dynamique
        while (tasks_completed < num_tasks) {

            std::cout << "Remaining tasks: " << num_tasks - tasks_completed << " ..." << std::endl;

            int task_index = -1;
            AggregatedResult result;
            double task_duration = 0.0;

            // Réception de l'indice de tâche correspondant
            MPI_Recv(&task_index, 1, MPI_INT, MPI_ANY_SOURCE, TAG_RESULT, MPI_COMM_WORLD, &status);
            int worker = status.MPI_SOURCE;

            // Réception d'un résultat et de la durée de la tâche
            MPI_Recv(&result, sizeof(AggregatedResult), MPI_BYTE, worker, TAG_RESULT, MPI_COMM_WORLD, &status);
            MPI_Recv(&task_duration, 1, MPI_DOUBLE, worker, TAG_RESULT, MPI_COMM_WORLD, &status);
            
            task_indices[task_index] = task_index;
            results[task_index] = result;
            durations[task_index].index = task_index;
            durations[task_index].rank = worker;
            durations[task_index].duration = task_duration;
 
            tasks_completed++;

            // S'il reste des tâches, on envoie la suivante à ce worker
            if (next_task < num_tasks) {
                MPI_Send(&next_task, 1, MPI_INT, worker, TAG_WORK, MPI_COMM_WORLD);
                next_task++;
            } else {
                // Sinon, on envoie un message d'arrêt
                MPI_Send(nullptr, 0, MPI_INT, worker, TAG_STOP, MPI_COMM_WORLD);
            }

            //
            // USLEEP ? nécessaire si je fais pas MPI_Recv mais MPI_Irecv ?
        }

        if (num_tasks < size - 1) { // -1 car le master n'est pas un worker
            // Si il y a plus de workers que de tâches, on doit envoyer un message d'arrêt aux workers inactifs
            for (int worker = num_tasks + 1; worker < size; ++worker) {
                MPI_Send(nullptr, 0, MPI_INT, worker, TAG_STOP, MPI_COMM_WORLD);
            }
        }

        auto end_time =  std::chrono::high_resolution_clock::now();
        double total_duration = std::chrono::duration<double> (end_time - start_time).count();

        std::cout << "All tasks completed !\nTotal duration: " << total_duration << " seconds." << std::endl;
        // Afficher les résultats

        std::ofstream result_file;
        std::string result_file_name = cfg.output.directory + "/aggregated_results.txt";
        result_file.open(result_file_name);
        result_file.precision(10);
        config::Config param_cfg = cfg; // copie de la config pour modifier les paramètres de chaque tâche

        result_file << "task" << "\t" << "#run" << "\t" << "N" << "\t" << "sat" << "\t" << "q" << "\t" 
                    << "mean_Nc" << "\t" << "mean_KE"  << "\t" << "duration" << std::endl;

        for (int i = 0; i < num_tasks; ++i) {
            // Récupérer les paramètres correspondants à la tâche i
            config::get_parameters(param_cfg, task_indices[i]);

            result_file << task_indices[i] << "\t"
                        << param_cfg.simulation.current_run << "\t"
                        << param_cfg.particles.N << "\t"
                        << param_cfg.patch.sat << "\t"
                        << param_cfg.patch.q << "\t"
                        << results[i].mean_caged_particles << "\t"
                        << results[i].mean_kinetic_energy << "\t"
                        << durations[i].duration << "\n";
        }
        result_file.close();

        std::vector<RankInfo> durations_by_rank(size-1, {-1, 0, 0.0, 0.0}); // size-1 car le master n'est pas un worker


        for (const auto& t : durations) {
            durations_by_rank[t.rank-1].total_duration += t.duration;
            durations_by_rank[t.rank-1].num_tasks += 1;
            durations_by_rank[t.rank-1].rank = t.rank;
            durations_by_rank[t.rank-1].idle_time = total_duration - durations_by_rank[t.rank-1].total_duration;
        }
    
        std::ofstream logs_file;
        std::string logs_file_name = cfg.output.directory + "/logs.txt";
        logs_file.open(logs_file_name);
        logs_file.precision(10);

        logs_file << "Rank" << "\t" << "Number of Tasks" << "\t" << "Duration (s)" << "\t" << "Idle time (s)" << std::endl;

        for (const auto& info : durations_by_rank) {
            logs_file << info.rank << "\t" << info.num_tasks << "\t\t" << info.total_duration << "\t" << info.idle_time << std::endl;
        }

        logs_file << "\n" << "Total duration: " << total_duration << " seconds." << std::endl;

        logs_file.close();

    } else {
        // ----- WORKER -----
        MPI_Status status;

        while (true) {
            int task_index = -1;
            MPI_Recv(&task_index, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if (status.MPI_TAG == TAG_STOP) 
                break; // sortir de la boucle

            config::get_parameters(cfg, task_index);

            // POUR LOG std::cout << "Process " << rank << " received task " << task_index << " param: " << "current_run=" << cfg.simulation.current_run << " N=" << cfg.particles.N << " sat=" << cfg.patch.sat << " q=" << cfg.patch.q << std::endl; 
            // std::cout << "Process " << rank << " received task " << task_index << " param: " << "current_run=" << cfg.simulation.current_run << " N=" << cfg.particles.N << " sat=" << cfg.patch.sat << " q=" << cfg.patch.q << std::endl; 

            auto task_start_time = std::chrono::high_resolution_clock::now();
            AggregatedResult res = simulation(cfg, task_index);
            auto task_end_time =  std::chrono::high_resolution_clock::now();
            double task_duration = std::chrono::duration<double> (task_end_time - task_start_time).count();
     

            // notifier le master que la tâche est terminée
            MPI_Send(&task_index, 1, MPI_INT, 0, TAG_RESULT, MPI_COMM_WORLD);
            MPI_Send(&res, sizeof(AggregatedResult), MPI_BYTE, 0, TAG_RESULT, MPI_COMM_WORLD);
            MPI_Send(&task_duration, 1, MPI_DOUBLE, 0, TAG_RESULT, MPI_COMM_WORLD);

            // std::cout << "Process " << rank << " sent completion notification for task " << task_index << std::endl;
            // POUR LOG : std::cout << "Process " << rank << " sent completion notification for task " << task_index << std::endl;
        }
    }


    MPI_Finalize(); // Fermeture de MPI
    return 0;
}