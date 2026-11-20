#include "reference/UniformFieldSolution.hpp"
#include "aurora/physics/Gyromotion.hpp"
#include "aurora/physics/Constants.hpp"
#include "aurora/physics/Relativity.hpp"
#include "aurora/simulation/TrajectoryRunner.hpp"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <optional>
#include <string>

namespace {
constexpr double tau = 2.0 * std::numbers::pi;

struct Measurement {
    double dt{}, position_error{}, momentum_error{}, energy_error{}, phase_error{}, period_error{},
        radius_error{};
};

Measurement measure(std::size_t per_period, std::size_t periods) {
    using namespace aurora;
    const Vec3 b{0, 0, 3e-5};
    const UniformField field(b);
    const ParticleState initial{
        {0, 0, 0}, momentum_from_energy(constants::eV_to_J(10000), {std::sqrt(3.0), 0, 1}), 0};
    const double period = gyroperiod_s(initial.momentum_kg_m_s, b);
    const double radius = larmor_radius_m(initial.momentum_kg_m_s, b);
    Measurement m;
    m.dt = period / static_cast<double>(per_period);
    double previous_phase = 0, accumulated_phase = 0;
    Vec3 quarter{}, half{};
    const auto summary = run_trajectory(initial, field, {m.dt, per_period * periods, 1},
        [&](std::size_t step, const ParticleState& s) {
            const auto exact = test_reference::uniform_field_solution(initial, b, s.time_s);
            m.position_error = std::max(m.position_error, norm(s.position_m - exact.position_m));
            m.momentum_error = std::max(m.momentum_error,
                norm(s.momentum_kg_m_s - exact.momentum_kg_m_s) / norm(initial.momentum_kg_m_s));
            const double phase = std::atan2(s.momentum_kg_m_s.y, s.momentum_kg_m_s.x);
            accumulated_phase += std::remainder(phase - previous_phase, tau);
            previous_phase = phase;
            if (step == per_period / 4)
                quarter = s.position_m;
            if (step == per_period / 2)
                half = s.position_m;
        });
    m.energy_error = summary.max_relative_energy_error;
    m.phase_error = std::abs(accumulated_phase - tau * static_cast<double>(periods));
    const double measured_period = tau * summary.final_state.time_s / accumulated_phase;
    m.period_error = std::abs(measured_period / period - 1.0);
    // Circumradius of three projected trajectory points, independent of the pusher formula.
    quarter.z = 0;
    half.z = 0;
    const double area_twice = norm(cross(quarter, half));
    const double measured_radius =
        norm(quarter) * norm(half) * norm(half - quarter) / (2.0 * area_twice);
    m.radius_error = std::abs(measured_radius / radius - 1.0);
    return m;
}
}

int main(int argc, char** argv) {
    try {
        std::filesystem::path output = "results/uniform_field_convergence.csv";
        if (argc == 3 && std::string(argv[1]) == "--output")
            output = argv[2];
        else if (argc != 1)
            throw std::invalid_argument(
                "Usage: aurora_uniform_field_convergence [--output results.csv]");

        if (output.has_parent_path())
            std::filesystem::create_directories(output.parent_path());
        std::ofstream csv(output);
        csv.exceptions(std::ios::failbit | std::ios::badbit);
        csv << std::setprecision(17)
            << "energy_keV,field_T,pitch_deg,steps_per_period,periods,dt_s,"
               "max_position_error_m,max_relative_momentum_error,phase_error_rad,"
               "relative_period_error,relative_radius_error,max_relative_energy_error,"
               "position_order,momentum_order,passed\n";

        auto write_row = [&](std::size_t steps, std::size_t periods, const Measurement& m,
                             bool passed, std::optional<double> position_order = {},
                             std::optional<double> momentum_order = {}) {
            csv << 10 << ',' << 3e-5 << ',' << 60 << ',' << steps << ',' << periods << ',' << m.dt
                << ',' << m.position_error << ',' << m.momentum_error << ',' << m.phase_error << ','
                << m.period_error << ',' << m.radius_error << ',' << m.energy_error << ',';
            if (position_order)
                csv << *position_order;
            csv << ',';
            if (momentum_order)
                csv << *momentum_order;
            csv << ',' << (passed ? "true" : "false") << '\n';
        };

        bool passed = true;
        Measurement previous{};
        for (const std::size_t n : {40u, 80u, 160u, 320u}) {
            const auto m = measure(n, 1);
            bool row_passed = std::isfinite(m.radius_error) && m.energy_error < 1e-9;
            std::optional<double> position_order, momentum_order;
            if (previous.position_error > 0) {
                position_order = std::log2(previous.position_error / m.position_error);
                momentum_order = std::log2(previous.momentum_error / m.momentum_error);
                row_passed = row_passed && *position_order > 1.9 && *position_order < 2.1 &&
                    *momentum_order > 1.9 && *momentum_order < 2.1;
            }
            write_row(n, 1, m, row_passed, position_order, momentum_order);
            passed = passed && row_passed;
            previous = m;
        }

        const auto long_run = measure(200, 100);
        const bool long_run_passed = long_run.energy_error < 1e-9 && long_run.phase_error > 0 &&
            long_run.phase_error < 0.06 && long_run.period_error < 1e-4 &&
            long_run.radius_error < 2e-4;
        write_row(200, 100, long_run, long_run_passed);
        passed = passed && long_run_passed;
        csv.close();

        std::cout << (passed ? "PASS: " : "FAIL: ") << output.string() << '\n';
        return passed ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "aurora_uniform_field_convergence: " << e.what() << '\n';
        return 1;
    }
}
