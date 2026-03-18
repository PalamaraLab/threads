#ifndef THREADS_ARG_THREADS_IO_HPP
#define THREADS_ARG_THREADS_IO_HPP

#include "ThreadingInstructions.hpp"
#include <string>
#include <vector>

void serialize_threads(
    const std::string& filename,
    ThreadingInstructions& instructions,
    const std::vector<std::vector<std::string>>& metadata_cols,
    const std::vector<double>& allele_ages,
    const std::vector<std::string>& sample_names);

ThreadingInstructions deserialize_threads(const std::string& filename);

// Read optional string datasets
std::vector<std::vector<std::string>> read_threads_metadata(const std::string& filename);
std::vector<std::string> read_threads_sample_names(const std::string& filename);

#endif
