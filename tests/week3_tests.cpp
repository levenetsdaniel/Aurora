#include "aurora/io/ExperimentFiles.hpp"
#include "aurora/physics/Relativity.hpp"
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>

namespace {
using namespace aurora;

void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

template <class F> void rejects(F f) {
    try {
        f();
    } catch (const std::exception&) {
        return;
    }
    throw std::runtime_error("Expected rejection");
}

void classification() {
    RunSummary r;
    const std::pair<StopReason, ExperimentOutcome> cases[] = {
        {StopReason::InnerBoundary, ExperimentOutcome::AtmosphereReached},
        {StopReason::Mirror, ExperimentOutcome::Reflected},
        {StopReason::OuterBoundary, ExperimentOutcome::Escaped},
        {StopReason::TimeLimit, ExperimentOutcome::TimeLimit},
        {StopReason::StepLimit, ExperimentOutcome::StepLimit},
        {StopReason::NumericalError, ExperimentOutcome::NumericalError}};
    for (auto [reason, expected] : cases) {
        r.stop_reason = reason;
        require(classify_atmosphere_outcome(r) == expected, "Wrong outcome classification");
    }
    r.events.push_back({EventType::Mirror, {}, 0});
    r.stop_reason = StopReason::TimeLimit;
    require(classify_atmosphere_outcome(r) == ExperimentOutcome::TimeLimit,
        "A limited run with mirror history must not be called a completed reflection experiment");
    for (auto o : {ExperimentOutcome::TimeLimit, ExperimentOutcome::StepLimit,
             ExperimentOutcome::NumericalError}) {
        require(!is_resolved(o), "Unfinished/error outcome was counted as resolved");
        require(compare_prediction(LossConePrediction::Reflection, o) ==
                PredictionAgreement::Unresolved,
            "Unresolved calculation was treated as a prediction mismatch or match");
    }
    require(compare_prediction(LossConePrediction::AtmosphereEntry,
                ExperimentOutcome::AtmosphereReached) == PredictionAgreement::Match,
        "Entry match");
    require(compare_prediction(LossConePrediction::Reflection,
                ExperimentOutcome::AtmosphereReached) == PredictionAgreement::Mismatch,
        "Entry mismatch");
    require(compare_prediction(LossConePrediction::Boundary, ExperimentOutcome::Reflected) ==
            PredictionAgreement::Boundary,
        "Boundary must be separate");
    require(compare_prediction(std::nullopt, ExperimentOutcome::Reflected) ==
            PredictionAgreement::NotApplicable,
        "Absent prediction");
    require(compare_prediction(LossConePrediction::Reflection, ExperimentOutcome::Escaped) ==
            PredictionAgreement::NotApplicable,
        "Escaped comparison");
}

void parameters() {
    ExperimentParameters p;
    p.energy_keV = 123.456789;
    p.pitch_deg = 175;
    p.gyro_phase_deg = 23;
    p.dipole_axis = {1, 2, -3};
    p.max_steps = 12345678;
    p.event_tolerance_s = 3e-12;
    p.stop_at_first_mirror = false;
    std::ostringstream first;
    write_experiment_parameters(first, p);
    std::istringstream input(first.str());
    const auto loaded = read_experiment_parameters(input);
    std::ostringstream second;
    write_experiment_parameters(second, loaded);
    require(first.str() == second.str(), "Parameter round trip lost values or precision");
    const auto prepared = prepare_experiment(p);
    require(std::abs(dot(normalized(prepared.initial.position_m), prepared.field.axis())) < 1e-14,
        "Tilted equatorial start");
    require(std::abs(kinetic_energy_J(prepared.initial.momentum_kg_m_s) /
                    constants::eV_to_J(p.energy_keV * 1000) -
                1) < 1e-14,
        "Initial energy");
    for (const std::string suffix : {"\nunknown=1\n", "\nenergy_keV=1\n"}) {
        std::istringstream bad(first.str() + suffix);
        rejects([&] {
            read_experiment_parameters(bad);
        });
    }
    std::istringstream missing("schema_version=1\n");
    rejects([&] {
        read_experiment_parameters(missing);
    });
    for (double h : {-1., 0., std::numeric_limits<double>::quiet_NaN()}) {
        auto bad = p;
        bad.atmosphere_height_m = h;
        rejects([&] {
            prepare_experiment(bad);
        });
    }
    for (double pitch : {-1., 181., std::numeric_limits<double>::infinity()}) {
        auto bad = p;
        bad.pitch_deg = pitch;
        rejects([&] {
            prepare_experiment(bad);
        });
    }
    auto bad = p;
    bad.equatorial_distance_RE = 1.001;
    rejects([&] {
        prepare_experiment(bad);
    });
    bad = p;
    bad.outer_distance_RE = 1;
    rejects([&] {
        prepare_experiment(bad);
    });
    bad = p;
    bad.dt_s = 0;
    rejects([&] {
        prepare_experiment(bad);
    });
    bad = p;
    bad.dipole_axis = {};
    rejects([&] {
        prepare_experiment(bad);
    });
}

void hemispheres() {
    for (double pitch : {5., 45., 135., 175.}) {
        ExperimentParameters p;
        p.pitch_deg = pitch;
        p.dt_s = 4e-8;
        const auto e = prepare_experiment(p);
        const auto r = run_atmosphere_experiment(e);
        const auto expected = (pitch == 5 || pitch == 175) ? ExperimentOutcome::AtmosphereReached
                                                           : ExperimentOutcome::Reflected;
        require(r.outcome == expected && r.agreement == PredictionAgreement::Match,
            "Both hemisphere prediction check");
        require(
            r.trajectory.final_state.position_m.z * (pitch < 90 ? 1 : -1) > 0, "Wrong hemisphere");
        require(r.trajectory.max_relative_energy_error < 1e-9, "Energy drift");
    }
    ExperimentParameters p;
    p.max_time_s = 1e-6;
    require(
        run_atmosphere_experiment(prepare_experiment(p)).outcome == ExperimentOutcome::TimeLimit,
        "Time limit");
    p.max_time_s = 0.1;
    p.max_steps = 0;
    require(
        run_atmosphere_experiment(prepare_experiment(p)).outcome == ExperimentOutcome::StepLimit,
        "Zero steps");
    p.max_steps = 100;
    p.outer_distance_RE = 2.000001;
    require(run_atmosphere_experiment(prepare_experiment(p)).outcome == ExperimentOutcome::Escaped,
        "Outer escape");
    p = ExperimentParameters{};
    p.energy_keV = 0;
    p.max_time_s = 1e-6;
    const auto stationary = prepare_experiment(p);
    require(!stationary.prediction, "Rest particle has no directional prediction");
    require(run_atmosphere_experiment(stationary).outcome == ExperimentOutcome::TimeLimit,
        "Rest stays unresolved");
}

void reproducibility() {
    const auto root = std::filesystem::path("test-results/reproducibility");
    ExperimentParameters p;
    p.max_steps = 10;
    p.sample_every = 3;
    const auto first = save_experiment(p, root);
    std::ifstream parameters_file(first.directory / "parameters.txt");
    const auto loaded = read_experiment_parameters(parameters_file);
    const auto second = save_experiment(loaded, root);
    require(first.directory != second.directory, "Results were overwritten");
    require(first.result.outcome == ExperimentOutcome::StepLimit &&
            second.result.agreement == PredictionAgreement::Unresolved,
        "Incomplete result status");
    require(norm(first.result.trajectory.final_state.position_m -
                second.result.trajectory.final_state.position_m) == 0,
        "Replay position differs");
    require(norm(first.result.trajectory.final_state.momentum_kg_m_s -
                second.result.trajectory.final_state.momentum_kg_m_s) == 0,
        "Replay momentum differs");
    for (const char* name : {"parameters.txt", "manifest.txt", "summary.csv", "events.csv",
             "trajectory.csv", "source/source.sha256", "source/main.cpp", "source/CMakeLists.txt"})
        require(std::filesystem::is_regular_file(first.directory / name),
            "Missing reproducibility artifact");
    std::ifstream manifest(first.directory / "manifest.txt");
    const std::string text((std::istreambuf_iterator<char>(manifest)), {});
    require(text.find("execution_status=finished") != std::string::npos &&
            text.find("source_sha256=") != std::string::npos,
        "Manifest status/provenance");
    auto invalid = p;
    invalid.atmosphere_height_m = -1;
    const auto invalid_root = first.directory / "invalid";
    rejects([&] {
        save_experiment(invalid, invalid_root);
    });
    require(!std::filesystem::exists(invalid_root), "Invalid configuration created output files");
    std::ostringstream failed;
    failed.setstate(std::ios::badbit);
    rejects([&] {
        write_experiment_parameters(failed, p);
    });
}
}

int main(int argc, char** argv) {
    const std::map<std::string, void (*)()> suites{{"classification", classification},
        {"parameters", parameters}, {"hemispheres", hemispheres},
        {"reproducibility", reproducibility}};
    try {
        if (argc != 2 || !suites.contains(argv[1]))
            throw std::invalid_argument("Unknown week-three test suite");
        suites.at(argv[1])();
        std::cout << "PASS: " << argv[1] << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FAIL: " << e.what() << '\n';
        return 1;
    }
}
