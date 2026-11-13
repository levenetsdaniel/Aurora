#pragma once

#include "aurora/physics/LossCone.hpp"
#include "aurora/simulation/TrajectoryRunner.hpp"

namespace aurora {
struct ExperimentParameters {
    double energy_keV{100};
    double pitch_deg{45};
    double gyro_phase_deg{0};
    double equatorial_distance_RE{2};
    double dt_s{2e-8};
    double max_time_s{0.1};
    double atmosphere_height_m{100000};
    double outer_distance_RE{4};
    double dipole_equatorial_T{3.12e-5};
    Vec3 dipole_axis{0, 0, -1};
    std::size_t max_steps{10000000};
    std::size_t sample_every{100};
    double event_tolerance_s{1e-11};
    double mirror_hysteresis{1e-3};
    bool stop_at_first_mirror{true};
};

enum class ExperimentOutcome {
    AtmosphereReached,
    Reflected,
    Escaped,
    TimeLimit,
    StepLimit,
    NumericalError
};
enum class PredictionAgreement { Match, Mismatch, Boundary, Unresolved, NotApplicable };
const char* outcome_name(ExperimentOutcome outcome);
const char* agreement_name(PredictionAgreement agreement);
const char* prediction_name(LossConePrediction prediction);
bool is_resolved(ExperimentOutcome outcome);
// Inner boundary must represent the top of the atmosphere in this experiment.
ExperimentOutcome classify_atmosphere_outcome(const RunSummary& run);
PredictionAgreement compare_prediction(
    std::optional<LossConePrediction> prediction, ExperimentOutcome outcome);

struct PreparedExperiment {
    DipoleField field;
    ParticleState initial;
    RunConfig config;
    EquatorialLossCone loss_cone;
    std::optional<LossConePrediction> prediction;
};

struct ExperimentResult {
    RunSummary trajectory;
    ExperimentOutcome outcome;
    PredictionAgreement agreement;
};

// Validates every parameter before the caller creates output files.
PreparedExperiment prepare_experiment(const ExperimentParameters& parameters);
ExperimentResult run_atmosphere_experiment(
    const PreparedExperiment& experiment, const TrajectoryObserver& observer = {});
} // namespace aurora
