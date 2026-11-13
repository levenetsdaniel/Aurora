#include "aurora/fields/DipoleField.hpp"
#include "aurora/physics/ParticleProperties.hpp"
#include "aurora/physics/Gyromotion.hpp"
#include "aurora/simulation/TrajectoryRunner.hpp"
#include "aurora/simulation/EventLocator.hpp"
#include "aurora/io/TrajectoryCsv.hpp"
#include "aurora/physics/Relativity.hpp"
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <numbers>
#include <sstream>
#include <vector>

namespace {
using namespace aurora;
constexpr double R = constants::earth_radius_m;

void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance, const char* label) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
        std::ostringstream out;
        out.precision(17);
        out << label << ": " << actual << " vs " << expected << " +/- " << tolerance;
        throw std::runtime_error(out.str());
    }
}

template <class F> void rejects(F f) {
    try {
        f();
    } catch (const std::exception&) {
        return;
    }
    throw std::runtime_error("Expected invalid input to be rejected");
}

ParticleState electron(double keV = 100, double pitch_deg = 45) {
    const double a = pitch_deg * std::numbers::pi / 180;
    return {{2 * R, 0, 0},
        momentum_from_energy(constants::eV_to_J(keV * 1000), {std::sin(a), 0, std::cos(a)}), 0};
}

void dipole() {
    DipoleField f;
    near(norm(f.at({R, 0, 0}, 0) - Vec3{0, 0, 3.12e-5}), 0, 1e-19, "equatorial field");
    near(norm(f.at({0, 0, R}, 0) - Vec3{0, 0, -6.24e-5}), 0, 1e-19, "polar field");
    const Vec3 p{1.3 * R, 0.7 * R, -0.8 * R};
    near(norm(f.at(p * 2, 0)) / norm(f.at(p, 0)), 0.125, 1e-16, "inverse cube");
    near(norm(f.at(-p, 0) - f.at(p, 0)), 0, 1e-19, "inversion symmetry");
    const DipoleField tilted(3.12e-5, {1, 0, 0});
    near(norm(tilted.at({0, R, 0}, 0) - Vec3{-3.12e-5, 0, 0}), 0, 1e-19, "tilted equator");
    near(norm(tilted.at({R, 0, 0}, 0) - Vec3{6.24e-5, 0, 0}), 0, 1e-19, "tilted pole");
    const double h = R * 1e-5;
    const double divergence = (f.at(p + Vec3{h, 0, 0}, 0).x - f.at(p - Vec3{h, 0, 0}, 0).x +
                                  f.at(p + Vec3{0, h, 0}, 0).y - f.at(p - Vec3{0, h, 0}, 0).y +
                                  f.at(p + Vec3{0, 0, h}, 0).z - f.at(p - Vec3{0, 0, h}, 0).z) /
        (2 * h);
    near(divergence / (norm(f.at(p, 0)) / R), 0, 1e-8, "div B");
    Vec3 gradient;
    gradient.x = (norm(f.at(p + Vec3{h, 0, 0}, 0)) - norm(f.at(p - Vec3{h, 0, 0}, 0))) / (2 * h);
    gradient.y = (norm(f.at(p + Vec3{0, h, 0}, 0)) - norm(f.at(p - Vec3{0, h, 0}, 0))) / (2 * h);
    gradient.z = (norm(f.at(p + Vec3{0, 0, h}, 0)) - norm(f.at(p - Vec3{0, 0, h}, 0))) / (2 * h);
    near(norm(gradient - f.magnitude_gradient(p)) / norm(gradient), 0, 1e-8, "gradient B");
    for (double theta : {0.8, 1.1, 1.57, 2.0, 2.3}) {
        const double r = 3 * R * std::pow(std::sin(theta), 2);
        const double dr = 6 * R * std::sin(theta) * std::cos(theta);
        const Vec3 position{r * std::sin(theta), 0, r * std::cos(theta)};
        const Vec3 tangent{dr * std::sin(theta) + r * std::cos(theta), 0,
            dr * std::cos(theta) - r * std::sin(theta)};
        near(norm(cross(normalized(tangent), normalized(f.at(position, 0)))), 0, 1e-14,
            "field line r=L R sin^2(theta)");
    }
    rejects([&] {
        f.at({}, 0);
    });
    rejects([&] {
        f.at({0.49 * R, 0, 0}, 0);
    });
    rejects([] {
        DipoleField invalid(-1);
    });
    rejects([] {
        DipoleField invalid(1, {});
    });
    rejects([] {
        DipoleField invalid(1, {0, 0, 1}, R, 2 * R);
    });
}

void particle_properties_test() {
    const UniformField b({0, 0, 2e-5});
    auto s = electron(100, 30);
    auto d = particle_properties(s, b);
    near(*d.pitch_angle_rad, std::numbers::pi / 6, 1e-15, "pitch angle");
    const double p = norm(s.momentum_kg_m_s);
    near(d.parallel_momentum_kg_m_s / p, std::sqrt(3.0) / 2, 1e-15, "parallel momentum");
    near(
        d.magnetic_moment_J_T / (p * p / (8 * constants::electron_mass_kg * 2e-5)), 1, 1e-14, "mu");
    near(predicted_mirror_field_T(s, b), 8e-5, 1e-18, "B mirror");
    s.momentum_kg_m_s = -s.momentum_kg_m_s;
    d = particle_properties(s, b);
    near(*d.pitch_angle_rad, 5 * std::numbers::pi / 6, 1e-15, "backwards pitch");
    s.momentum_kg_m_s = {};
    require(!particle_properties(s, b).pitch_angle_rad, "pitch undefined at rest");
    near(particle_properties(s, b).magnetic_moment_J_T, 0, 0, "mu at rest");
    rejects([&] {
        particle_properties(s, ZeroField{});
    });
    rejects([&] {
        predicted_mirror_field_T(s, b);
    });
    std::ostringstream csv;
    write_magnetic_trajectory_header(csv);
    write_magnetic_trajectory_row(csv, 0, electron(), DipoleField{});
    require(csv.str().find("pitch_deg,p_parallel_kg_m_s,mu_J_T") != std::string::npos,
        "magnetic CSV columns");
}

void events() {
    const ZeroField field;
    const ParticleState inward{{2, 0, 0}, momentum_from_velocity({-1000, 0, 0}), 0};
    for (double dt : {0.004, 0.003, 0.0017}) {
        RunConfig c{dt, 10, 10};
        c.inner_radius_m = 1;
        c.outer_radius_m = 5;
        c.event_tolerance_s = 1e-13;
        const auto s = run_trajectory(inward, field, c);
        require(s.stop_reason == StopReason::InnerBoundary, "inward crossing status");
        near(s.final_state.time_s, 0.001, 1e-13,
            "FIRST inward crossing, even enter-exit in one step");
        near(s.final_state.position_m.x, 1, 1e-10, "inward crossing position");
        require(s.events.size() == 1 && s.events[0].bracket_width_s <= c.event_tolerance_s,
            "localized boundary");
    }
    const ParticleState outward{{2, 0, 0}, momentum_from_velocity({1000, 0, 0}), 0};
    RunConfig c{0.004, 10, 10};
    c.inner_radius_m = 1;
    c.outer_radius_m = 3;
    const auto out = run_trajectory(outward, field, c);
    require(out.stop_reason == StopReason::OuterBoundary, "outward status");
    near(out.final_state.time_s, 0.001, 1e-11, "outward time");
    // A curved trajectory: independently known intersection for x=rL*sin(omega*t).
    const UniformField uniform({0, 0, 3e-5});
    const ParticleState start{{}, momentum_from_energy(constants::eV_to_J(10000), {1, 0, 0}), 0};
    const double period = gyroperiod_s(start.momentum_kg_m_s, uniform.value_T());
    const double radius = larmor_radius_m(start.momentum_kg_m_s, uniform.value_T());
    double previous_error = 0;
    for (std::size_t steps : {40u, 80u, 160u}) {
        RunConfig config{period / static_cast<double>(steps), steps, 7};
        config.outer_radius_m = radius;
        config.event_tolerance_s = period * 1e-10;
        const auto s = run_trajectory(start, uniform, config);
        require(s.stop_reason == StopReason::OuterBoundary, "curved boundary status");
        const double error = std::abs(s.final_state.time_s - period / 6);
        if (previous_error > 0)
            require(previous_error / error > 3 && previous_error / error < 5,
                "event time second-order convergence");
        previous_error = error;
    }
    // Overshooting the dipole cutoff must not evaluate the field there.
    const DipoleField dipole;
    const ParticleState falling{{1.1 * R, 0, 0}, momentum_from_velocity({-1000, 0, 0}), 0};
    RunConfig coarse{10000, 1, 1};
    coarse.inner_radius_m = R;

    // Use a field with an explicit inner domain, but zero force, to isolate the boundary algorithm.
    class CutoffField final : public MagneticField {
        Vec3 at(Vec3 p, double) const override {
            if (norm(p) < 0.5 * R)
                throw std::domain_error("cutoff");
            return {};
        }
    } cutoff;

    const auto crossing = run_trajectory(falling, cutoff, coarse);
    require(crossing.stop_reason == StopReason::InnerBoundary,
        "localize before endpoint field evaluation");
    near(norm(crossing.final_state.position_m) / R, 1, 1e-12, "cutoff boundary");
}

void limits() {
    const ZeroField field;
    const ParticleState s{{2, 0, 0}, momentum_from_velocity({1000, 0, 0}), 3};
    RunConfig c{0.01, 10, 4};
    c.max_duration_s = 0.025;
    std::vector<double> times;
    const auto result = run_trajectory(s, field, c, [&](std::size_t, const ParticleState& state) {
        times.push_back(state.time_s);
    });
    require(result.stop_reason == StopReason::TimeLimit, "time limit status");
    near(result.final_state.time_s, 3.025, 1e-15, "fractional final time");
    near(result.final_state.position_m.x, 27, 1e-10, "fractional final position");
    require(times.size() == 2 && times.back() == result.final_state.time_s,
        "terminal event exported despite stride");
    c.max_duration_s = 0;
    require(run_trajectory(s, field, c).completed_steps == 0, "zero duration");
    c = RunConfig{0.01, 2, 1};
    require(run_trajectory(s, field, c).stop_reason == StopReason::StepLimit, "step limit");
    c.inner_radius_m = 3;
    rejects([&] {
        run_trajectory(s, field, c);
    });
    c = RunConfig{0.01, 2, 1};
    c.outer_radius_m = -1;
    rejects([&] {
        run_trajectory(s, field, c);
    });
    c = RunConfig{0.01, 2, 1};
    c.stop_at_first_mirror = true;
    rejects([&] {
        run_trajectory(s, field, c);
    });
    c = RunConfig{0.01, 2, 1};
    c.event_tolerance_s = 0;
    rejects([&] {
        run_trajectory(s, field, c);
    });

    class FailingField final : public MagneticField {
        Vec3 at(Vec3, double t) const override {
            if (t > 3)
                throw std::runtime_error("field unavailable");
            return {};
        }
    } failing;

    c = RunConfig{0.01, 2, 1};
    const auto failed = run_trajectory(s, failing, c);
    require(failed.stop_reason == StopReason::NumericalError && !failed.error_message.empty(),
        "numerical failure status");
    near(failed.final_state.time_s, 3, 0, "failure returns last valid state");
    rejects([&] {
        run_trajectory(s, field, c, [](std::size_t, const ParticleState&) {
            throw std::runtime_error("I/O");
        });
    });
}

void dipole_convergence() {
    const DipoleField field;
    const auto initial = electron();
    const double duration = 0.001;
    const auto reference = run_trajectory(initial, field, {2.5e-9, 400000, 1}).final_state;
    double previous_r = 0, previous_p = 0;
    for (std::size_t steps : {25000u, 50000u, 100000u}) {
        const auto s =
            run_trajectory(initial, field, {duration / static_cast<double>(steps), steps, 1});
        const double error_r = norm(s.final_state.position_m - reference.position_m);
        const double error_p = norm(s.final_state.momentum_kg_m_s - reference.momentum_kg_m_s);
        if (previous_r > 0) {
            require(previous_r / error_r > 3.5 && previous_r / error_r < 5,
                "dipole position convergence");
            require(previous_p / error_p > 3.5 && previous_p / error_p < 5,
                "dipole momentum convergence");
        }
        previous_r = error_r;
        previous_p = error_p;
    }
}

void mirrors() {
    const DipoleField field;
    const auto initial = electron();
    RunConfig c{4e-8, 5000000, 1000};
    c.detect_mirrors = true;
    c.stop_at_first_mirror = true;
    c.inner_radius_m = 1.02 * R;
    c.outer_radius_m = 4 * R;
    const auto s = run_trajectory(initial, field, c);
    require(s.stop_reason == StopReason::Mirror && s.events.size() == 1, "first mirror");
    const double measured = norm(field.at(s.final_state.position_m, s.final_state.time_s));
    near(measured / predicted_mirror_field_T(initial, field), 1, 0.01,
        "adiabatic mirror prediction");
    near(particle_properties(s.final_state, field).parallel_fraction, 0, 1e-7,
        "zero parallel momentum");
    require(s.max_relative_magnetic_moment_error < 0.002, "first invariant in adiabatic dipole");
    require(s.max_relative_energy_error < 1e-9, "energy in dipole");
    c.stop_at_first_mirror = false;
    c.max_duration_s = 0.23;
    c.dt_s = 5e-8;
    const auto bounce = run_trajectory(initial, field, c);
    require(bounce.stop_reason == StopReason::TimeLimit && bounce.events.size() == 2,
        "two mirrors without duplicate zero crossings");
    require(bounce.events[0].state.position_m.z > 0 && bounce.events[1].state.position_m.z < 0,
        "opposite hemisphere mirrors");
    require(bounce.events[0].state.time_s < bounce.events[1].state.time_s, "ordered mirrors");
    const UniformField uniform({0, 0, 3e-5});
    c = RunConfig{1e-9, 100, 1};
    c.detect_mirrors = true;
    require(run_trajectory(
                {{}, momentum_from_energy(constants::eV_to_J(1e4), {1, 0, 0}), 0}, uniform, c)
                .events.empty(),
        "90-degree uniform orbit is not repeated reflection");
}
}

int main(int argc, char** argv) {
    const std::map<std::string, std::function<void()>> suites{{"dipole", dipole},
        {"particle_properties", particle_properties_test}, {"events", events}, {"limits", limits},
        {"mirrors", mirrors}, {"dipole_convergence", dipole_convergence}};
    try {
        if (argc != 2 || !suites.contains(argv[1]))
            throw std::invalid_argument("Unknown suite");
        suites.at(argv[1])();
        std::cout << "PASS: " << argv[1] << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
