#include "aurora/fields/DipoleField.hpp"
#include "aurora/io/EventCsv.hpp"
#include "aurora/io/TrajectoryCsv.hpp"
#include "aurora/physics/Relativity.hpp"
#include "aurora/physics/LossCone.hpp"
#include "aurora/simulation/TrajectoryRunner.hpp"
#include <iostream>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

namespace {
using namespace aurora;

void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

std::vector<std::string> columns(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream input(line);
    std::string value;
    while (std::getline(input, value, ','))
        result.push_back(value);
    return result;
}

void check_scenario(double pitch_deg, StopReason expected_reason, const char* expected_type) {
    const DipoleField field;
    const double pitch = pitch_deg * std::numbers::pi / 180;
    const ParticleState initial{{2 * constants::earth_radius_m, 0, 0},
        momentum_from_energy(constants::eV_to_J(100000), {std::sin(pitch), 0, std::cos(pitch)}), 0};
    RunConfig config{2e-8, 10000000, 1000};
    config.inner_radius_m = constants::earth_radius_m + 100000;
    config.outer_radius_m = 4 * constants::earth_radius_m;
    config.max_duration_s = 0.1;
    config.detect_mirrors = true;
    config.stop_at_first_mirror = true;
    const auto cone = equatorial_loss_cone(field, initial.position_m, *config.inner_radius_m);
    const auto predicted = predict_equatorial_outcome(cone, pitch);
    require(predicted ==
            (expected_reason == StopReason::InnerBoundary ? LossConePrediction::AtmosphereEntry
                                                          : LossConePrediction::Reflection),
        "Analytical loss cone prediction disagrees with the expected scenario");
    ParticleState last_exported{};
    std::size_t last_step = 0;
    std::ostringstream trajectory;
    write_magnetic_trajectory_header(trajectory);
    const auto run =
        run_trajectory(initial, field, config, [&](std::size_t step, const ParticleState& state) {
            last_exported = state;
            last_step = step;
            write_magnetic_trajectory_row(trajectory, step, state, field);
        });
    require(run.stop_reason == expected_reason, "Unexpected trajectory outcome");
    require(run.events.size() == 1, "Expected exactly one terminating event");
    const auto& event = run.events.front();
    require(event.state.time_s == run.final_state.time_s, "Final time differs from event time");
    require(norm(event.state.position_m - last_exported.position_m) == 0,
        "Trajectory export missed the localized event position");
    require(norm(event.state.momentum_kg_m_s - last_exported.momentum_kg_m_s) == 0,
        "Trajectory and event momenta differ");
    require(last_step == run.completed_steps, "Final step was not exported");
    require(run.max_relative_energy_error < 1e-9, "Energy is not conserved");

    std::ostringstream csv;
    write_event_header(csv);
    write_event_row(csv, event, field, "atmosphere_entry");
    std::istringstream input(csv.str());
    std::string header, row;
    std::getline(input, header);
    std::getline(input, row);
    const auto names = columns(header);
    const auto values = columns(row);
    require(names.size() == 12 && values.size() == names.size(), "Event CSV column count mismatch");
    require(values[0] == expected_type, "Wrong event name in CSV");
    require(names[7] == "altitude_m" && names[8] == "kinetic_energy_J" && names[9] == "px_kg_m_s" &&
            names[10] == "py_kg_m_s" && names[11] == "pz_kg_m_s",
        "Missing physical event columns");
    require(std::stod(values[1]) == last_exported.time_s, "CSV time lost precision");
    require(std::stod(values[8]) == kinetic_energy_J(event.state.momentum_kg_m_s),
        "Wrong event energy");
    require(std::stod(values[9]) == last_exported.momentum_kg_m_s.x &&
            std::stod(values[10]) == last_exported.momentum_kg_m_s.y &&
            std::stod(values[11]) == last_exported.momentum_kg_m_s.z,
        "Wrong exported momentum");
    if (expected_reason == StopReason::InnerBoundary) {
        require(event.type == EventType::InnerBoundary, "Entry must be an inner boundary event");
        require(std::abs(std::stod(values[7]) - 100000) < 0.005, "Entry altitude is not 100 km");
        require(event.bracket_width_s <= config.event_tolerance_s, "Entry was not localized");
        require(
            dot(event.state.position_m, velocity_from_momentum(event.state.momentum_kg_m_s)) < 0,
            "Atmosphere entry must be directed inward");
    } else {
        require(event.type == EventType::Mirror && std::stod(values[7]) > 100000,
            "Reflection must occur above the atmosphere");
    }
    std::cout << pitch_deg << " deg: " << expected_type << ", t=" << event.state.time_s
              << " s, altitude=" << values[7]
              << " m, energy=" << constants::J_to_eV(std::stod(values[8])) / 1000 << " keV\n";
}
}

int main() {
    try {
        check_scenario(45, StopReason::Mirror, "mirror");
        check_scenario(5, StopReason::InnerBoundary, "atmosphere_entry");
        std::ostringstream failed;
        failed.setstate(std::ios::badbit);
        bool rejected = false;
        try {
            aurora::write_event_header(failed);
        } catch (const std::exception&) {
            rejected = true;
        }
        require(rejected, "Event CSV write failure must be reported");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
