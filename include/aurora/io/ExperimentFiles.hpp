#pragma once

#include "aurora/simulation/AtmosphereExperiment.hpp"
#include <filesystem>
#include <iosfwd>

namespace aurora {
void write_experiment_parameters(std::ostream& out, const ExperimentParameters& parameters);
ExperimentParameters read_experiment_parameters(std::istream& input);

struct SavedExperiment {
    ExperimentResult result;
    std::filesystem::path directory;
    double wall_time_s;
};

// Creates a new directory for every run; existing results are never overwritten.
// The optional observer receives every step without changing the CSV sampling interval.
SavedExperiment save_experiment(const ExperimentParameters& parameters,
    const std::filesystem::path& output_root, const TrajectoryObserver& observer = {});
} // namespace aurora
