#pragma once

#include "aurora/numerics/RelativisticBorisPusher.hpp"
#include <functional>

namespace aurora {
struct RunConfig {
    double dt_s;
    std::size_t steps;
    std::size_t sample_every{1};
};
struct RunSummary {
    ParticleState final_state;
    std::size_t completed_steps{};
    double initial_energy_J{};
    double max_relative_energy_error{}; // 0 for an initially stationary electron.
    double max_speed_m_s{};
};
using TrajectoryObserver = std::function<void(std::size_t, const ParticleState&)>;
RunSummary run_trajectory(const ParticleState& initial, const MagneticField& field,
                          RunConfig config, const TrajectoryObserver& observer = {});
} // namespace aurora
