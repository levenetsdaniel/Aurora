#include "aurora/simulation/ElectronEnsemble.hpp"
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
using namespace aurora;

void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

template <class Action> void rejects(Action action) {
    try {
        action();
    } catch (const std::invalid_argument&) {
        return;
    }
    throw std::runtime_error("Invalid ensemble was accepted");
}

void ordering_and_single_particle_equivalence() {
    const std::vector<ElectronSample> samples{
        {23, 10, 15, 0}, {0, 100, 45, 90}, {88, 1000, 175, 270}};
    EnsembleParameters parameters;
    parameters.particle_count = samples.size();
    parameters.particle.max_steps = 4;

    const auto results = run_ensemble(parameters, samples);
    require(results.size() == samples.size(), "Wrong result count");
    for (std::size_t i = 0; i < samples.size(); ++i) {
        const auto& sample = samples[i];
        const auto& result = results[i];
        require(result.initial.id == sample.id && result.initial.energy_keV == sample.energy_keV &&
                result.initial.pitch_deg == sample.pitch_deg &&
                result.initial.gyro_phase_deg == sample.gyro_phase_deg,
            "Initial conditions or sample order changed");

        auto particle = parameters.particle;
        particle.energy_keV = sample.energy_keV;
        particle.pitch_deg = sample.pitch_deg;
        particle.gyro_phase_deg = sample.gyro_phase_deg;
        const auto single = run_atmosphere_experiment(prepare_experiment(particle));
        const auto& actual = result.simulation.trajectory;
        const auto& expected = single.trajectory;

        require(result.simulation.outcome == ExperimentOutcome::StepLimit &&
                result.simulation.agreement == single.agreement &&
                actual.stop_reason == expected.stop_reason &&
                actual.completed_steps == expected.completed_steps,
            "Ensemble changed the single-particle outcome");
        require(norm(actual.final_state.position_m - expected.final_state.position_m) == 0 &&
                norm(actual.final_state.momentum_kg_m_s - expected.final_state.momentum_kg_m_s) ==
                    0 &&
                actual.final_state.time_s == expected.final_state.time_s &&
                actual.max_relative_energy_error == expected.max_relative_energy_error &&
                actual.max_relative_magnetic_moment_error ==
                    expected.max_relative_magnetic_moment_error,
            "Ensemble calculation differs from a separate calculation");
    }
}

void known_physical_outcomes() {
    const std::vector<ElectronSample> samples{
        {42, 100, 5, 0}, {7, 100, 45, 90}, {100, 100, 135, 180}, {3, 100, 175, 270}};
    EnsembleParameters parameters;
    parameters.particle_count = samples.size();
    const auto results = run_ensemble(parameters, samples);

    for (std::size_t i = 0; i < results.size(); ++i) {
        const bool entry = i == 0 || i == 3;
        const auto& result = results[i].simulation;
        require(result.outcome ==
                (entry ? ExperimentOutcome::AtmosphereReached : ExperimentOutcome::Reflected),
            "Unexpected physical outcome");
        require(
            result.agreement == PredictionAgreement::Match && result.trajectory.events.size() == 1,
            "Missing physical event or analytical agreement");
        require(result.trajectory.events.front().type ==
                (entry ? EventType::InnerBoundary : EventType::Mirror),
            "Incorrect event type");
        require(result.trajectory.max_relative_energy_error < 1e-9, "Energy conservation failed");
    }
    require(results.front().simulation.trajectory.final_state.position_m.z > 0 &&
            results.back().simulation.trajectory.final_state.position_m.z < 0,
        "Entries must occur in opposite hemispheres");
}

void time_limit_and_continued_motion() {
    const std::vector<ElectronSample> samples{{5, 100, 45, 0}};
    EnsembleParameters parameters;
    parameters.particle_count = samples.size();
    parameters.particle.max_time_s = 0;
    auto results = run_ensemble(parameters, samples);
    require(results.front().simulation.outcome == ExperimentOutcome::TimeLimit &&
            results.front().simulation.trajectory.completed_steps == 0,
        "Time limit was relabeled as a resolved outcome");

    parameters.particle.max_time_s = 0.08;
    parameters.particle.stop_at_first_mirror = false;
    results = run_ensemble(parameters, samples);
    const auto& result = results.front().simulation;
    require(result.outcome == ExperimentOutcome::TimeLimit &&
            result.agreement == PredictionAgreement::Unresolved &&
            result.trajectory.events.size() == 1 &&
            result.trajectory.events.front().type == EventType::Mirror,
        "Ensemble ignored continued motion after reflection");
}

void invalid_inputs() {
    EnsembleParameters parameters;
    const std::vector<ElectronSample> valid{{1, 100, 45, 0}, {2, 100, 5, 0}};
    rejects([&] {
        run_ensemble(parameters, valid);
    });
    parameters.particle_count = 0;
    rejects([&] {
        run_ensemble(parameters, {});
    });
    parameters.particle_count = valid.size();

    auto invalid = valid;
    invalid.back().id = invalid.front().id;
    rejects([&] {
        run_ensemble(parameters, invalid);
    });
    for (double energy : {-1., std::numeric_limits<double>::quiet_NaN()}) {
        invalid = valid;
        invalid.back().energy_keV = energy;
        rejects([&] {
            run_ensemble(parameters, invalid);
        });
    }
    invalid = valid;
    invalid.back().gyro_phase_deg = std::numeric_limits<double>::infinity();
    rejects([&] {
        run_ensemble(parameters, invalid);
    });
    parameters.particle.dt_s = 0;
    rejects([&] {
        run_ensemble(parameters, valid);
    });

    // A huge step would fail during Boris initialization of the first particle.
    // Invalid input in the second particle must be detected before any integration starts.
    parameters.particle.dt_s = 1e200;
    parameters.particle.max_steps = 1;
    invalid = valid;
    invalid.back().pitch_deg = 181;
    rejects([&] {
        run_ensemble(parameters, invalid);
    });
}
}

int main() {
    try {
        ordering_and_single_particle_equivalence();
        known_physical_outcomes();
        time_limit_and_continued_motion();
        invalid_inputs();
        std::cout << "PASS: electron ensemble\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
