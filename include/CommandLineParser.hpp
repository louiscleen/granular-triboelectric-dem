#ifndef COMMANDLINEPARSER_HPP
#define COMMANDLINEPARSER_HPP

#include <string>

void parsing_options(int argc, char *argv[], std::string &config_file, std::string &seed,
                     std::string &output_path, int &exit_requested, int rank);

#endif
