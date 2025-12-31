#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <string_view>

inline constexpr std::string_view VERSION = "1.0.0";
inline constexpr std::string_view DEFAULT_CONFIG_FILE = "config.toml";

inline constexpr int TAG_WORK = 1;
inline constexpr int TAG_RESULT = 2;
inline constexpr int TAG_STOP = 3;

#endif
