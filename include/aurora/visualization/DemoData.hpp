#pragma once

#include "aurora/io/ExperimentFiles.hpp"
#include <vector>

namespace aurora::visualization {
struct VisualSample {
    double time_s{};
    Vec3 position_m;
    double altitude_km{};
    double parallel_fraction{};
};

struct DemoCase {
    ExperimentParameters parameters;
    std::filesystem::path directory;
    std::vector<VisualSample> trajectory;
    std::vector<VisualSample> gyro_segment;
    std::vector<TrajectoryEvent> events;
    ExperimentOutcome outcome;
    double energy_error{};
    double mu_error{};
};

DemoCase load_demo_case(const std::filesystem::path& directory);
std::vector<DemoCase> create_demo_cases(
    const ExperimentParameters& base, const std::filesystem::path& root);

// Display-only interpolation in physical time; it never advances the physics simulation.
VisualSample sample_at_time(const std::vector<VisualSample>& samples, double time_s);
std::vector<Vec3> dipole_field_line(const DipoleField& field, Vec3 equatorial_position_m,
    double atmosphere_radius_m, std::size_t points = 241);
} // namespace aurora::visualization
