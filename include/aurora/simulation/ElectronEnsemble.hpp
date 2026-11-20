#pragma once

#include "aurora/simulation/AtmosphereExperiment.hpp"
#include <cstddef>
#include <cstdint>
#include <vector>

namespace aurora {
struct ElectronSample {
    std::size_t id{};
    double energy_keV{};
    double pitch_deg{};
    double gyro_phase_deg{};
};

struct EnsembleParameters {
    ExperimentParameters particle;
    std::size_t particle_count{100};
    std::uint64_t seed{42};
};

struct ElectronResult {
    ElectronSample initial;
    ExperimentResult simulation;
};

std::vector<ElectronResult> run_ensemble(
    const EnsembleParameters& parameters, const std::vector<ElectronSample>& samples);
}
