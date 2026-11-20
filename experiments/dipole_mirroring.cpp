#include "aurora/fields/DipoleField.hpp"
#include "aurora/physics/ParticleProperties.hpp"
#include "aurora/physics/Relativity.hpp"
#include "aurora/simulation/TrajectoryRunner.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <optional>
#include <string_view>

namespace {
using namespace aurora;

struct Measurement {
    double energy_keV, pitch_deg, dt_s, seconds;
    ParticleState initial;
    RunSummary run;
    double mirror_error, latitude_deg, altitude_km;
};

Measurement measure(double keV, double pitch_deg, double dt) {
    const DipoleField field;
    const double pitch = pitch_deg * std::numbers::pi / 180;
    const ParticleState initial{{2 * constants::earth_radius_m, 0, 0},
        momentum_from_energy(constants::eV_to_J(keV * 1000), {std::sin(pitch), 0, std::cos(pitch)}),
        0};
    RunConfig c{dt, 100000000, 1000};
    c.max_duration_s = 1;
    c.detect_mirrors = true;
    c.stop_at_first_mirror = true;
    c.inner_radius_m = 1.02 * constants::earth_radius_m;
    c.outer_radius_m = 4 * constants::earth_radius_m;
    c.event_tolerance_s = 1e-13;
    const auto start = std::chrono::steady_clock::now();
    auto run = run_trajectory(initial, field, c);
    const double elapsed =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    if (run.stop_reason != StopReason::Mirror)
        throw std::runtime_error("No mirror: " + std::string(stop_reason_name(run.stop_reason)) +
            " " + run.error_message);
    const Vec3 p = run.final_state.position_m;
    const double error =
        std::abs(norm(field.at(p, 0)) / predicted_mirror_field_T(initial, field) - 1);
    return {keV, pitch_deg, dt, elapsed, initial, std::move(run), error,
        std::asin(p.z / norm(p)) * 180 / std::numbers::pi,
        (norm(p) - constants::earth_radius_m) / 1000};
}
}

int main(int argc, char** argv) {
    try {
        std::filesystem::path output = "results/dipole_mirroring.csv";
        if (argc == 3 && std::string(argv[1]) == "--output")
            output = argv[2];
        else if (argc != 1)
            throw std::invalid_argument("Usage: aurora_dipole_mirroring [--output results.csv]");

        if (output.has_parent_path())
            std::filesystem::create_directories(output.parent_path());
        std::ofstream csv(output);
        csv.exceptions(std::ios::failbit | std::ios::badbit);
        csv << std::setprecision(17)
            << "series,energy_keV,pitch_deg,dt_s,mirror_time_s,latitude_deg,altitude_km,"
               "measured_mirror_field_T,predicted_mirror_field_T,relative_mirror_field_error,"
               "max_relative_mu_error,max_relative_energy_error,reference_dt_s,"
               "time_difference_s,position_difference_m,event_bracket_width_s,steps,wall_time_s\n";

        const DipoleField field;
        auto write_row = [&](std::string_view series, const Measurement& m,
                             std::optional<double> reference_dt = {},
                             std::optional<double> time_error = {},
                             std::optional<double> position_error = {}) {
            csv << series << ',' << m.energy_keV << ',' << m.pitch_deg << ',' << m.dt_s << ','
                << m.run.final_state.time_s << ',' << m.latitude_deg << ',' << m.altitude_km << ','
                << norm(field.at(m.run.final_state.position_m, 0)) << ','
                << predicted_mirror_field_T(m.initial, field) << ',' << m.mirror_error << ','
                << m.run.max_relative_magnetic_moment_error << ','
                << m.run.max_relative_energy_error << ',';
            if (reference_dt)
                csv << *reference_dt;
            csv << ',';
            if (time_error)
                csv << *time_error;
            csv << ',';
            if (position_error)
                csv << *position_error;
            csv << ',' << m.run.events.front().bracket_width_s << ',' << m.run.completed_steps
                << ',' << m.seconds << '\n';
            csv.flush();
        };

        bool passed = true;
        for (double energy : {10., 100., 1000.})
            for (double angle : {30., 45., 60.}) {
                const auto m = measure(energy, angle, 4e-8);
                write_row("scan", m);
                passed = passed && m.mirror_error < 0.01 &&
                    m.run.max_relative_energy_error < 1e-9 &&
                    m.run.max_relative_magnetic_moment_error < 0.01;
            }

        const auto reference = measure(1000, 45, 2.5e-9);
        write_row("reference", reference, reference.dt_s, 0, 0);
        double previous_time = 0, previous_position = 0;
        for (double dt : {4e-8, 2e-8, 1e-8, 5e-9}) {
            const auto m = measure(1000, 45, dt);
            const double time_error =
                std::abs(m.run.final_state.time_s - reference.run.final_state.time_s);
            const double position_error =
                norm(m.run.final_state.position_m - reference.run.final_state.position_m);
            write_row("step_refinement", m, reference.dt_s, time_error, position_error);
            if (previous_time > 0)
                passed = passed && time_error < previous_time && position_error < previous_position;
            previous_time = time_error;
            previous_position = position_error;
        }
        csv.close();

        std::cout << (passed ? "PASS: " : "FAIL: ") << output.string() << '\n';
        return passed ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
