#pragma once

#include "aurora/numerics/RelativisticBorisPusher.hpp"
#include <functional>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace aurora {
struct RunConfig {
    double dt_s;
    std::size_t steps;
    std::size_t sample_every{1};
    std::optional<double> inner_radius_m{};
    std::optional<double> outer_radius_m{};
    double max_duration_s{std::numeric_limits<double>::infinity()};
    bool detect_mirrors{false};
    bool stop_at_first_mirror{false};
    double event_tolerance_s{1e-11};
    double mirror_hysteresis{1e-3}; // Fraction of |p|; rearm only after leaving the zero region.
};
enum class StopReason { StepLimit, TimeLimit, InnerBoundary, OuterBoundary, Mirror, NumericalError };
const char* stop_reason_name(StopReason reason);
enum class EventType { Mirror, InnerBoundary, OuterBoundary };
struct TrajectoryEvent {
    EventType type;
    ParticleState state;
    double bracket_width_s{};
};
struct RunSummary {
    ParticleState final_state;
    std::size_t completed_steps{};
    double initial_energy_J{};
    double max_relative_energy_error{}; // 0 for an initially stationary electron.
    double max_speed_m_s{};
    double initial_magnetic_moment_J_T{};
    double max_relative_magnetic_moment_error{};
    double max_field_T{};
    StopReason stop_reason{StopReason::StepLimit};
    std::vector<TrajectoryEvent> events{};
    std::string error_message{};
};
using TrajectoryObserver = std::function<void(std::size_t, const ParticleState&)>;
RunSummary run_trajectory(const ParticleState& initial, const MagneticField& field,
                          RunConfig config, const TrajectoryObserver& observer = {});
} // namespace aurora
