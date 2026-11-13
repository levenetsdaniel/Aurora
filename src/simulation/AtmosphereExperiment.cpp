#include "aurora/simulation/AtmosphereExperiment.hpp"
#include "aurora/physics/Relativity.hpp"
#include <numbers>

namespace aurora {
const char* outcome_name(ExperimentOutcome o) {
    switch (o) {
    case ExperimentOutcome::AtmosphereReached:
        return "atmosphere_reached";
    case ExperimentOutcome::Reflected:
        return "reflected";
    case ExperimentOutcome::Escaped:
        return "escaped";
    case ExperimentOutcome::TimeLimit:
        return "time_limit";
    case ExperimentOutcome::StepLimit:
        return "step_limit";
    case ExperimentOutcome::NumericalError:
        return "numerical_error";
    }
    return "unknown";
}

const char* agreement_name(PredictionAgreement a) {
    switch (a) {
    case PredictionAgreement::Match:
        return "match";
    case PredictionAgreement::Mismatch:
        return "mismatch";
    case PredictionAgreement::Boundary:
        return "boundary";
    case PredictionAgreement::Unresolved:
        return "unresolved";
    case PredictionAgreement::NotApplicable:
        return "not_applicable";
    }
    return "unknown";
}

const char* prediction_name(LossConePrediction p) {
    switch (p) {
    case LossConePrediction::AtmosphereEntry:
        return "atmosphere_entry";
    case LossConePrediction::Reflection:
        return "reflection";
    case LossConePrediction::Boundary:
        return "boundary";
    }
    return "unknown";
}

bool is_resolved(ExperimentOutcome o) {
    return o == ExperimentOutcome::AtmosphereReached || o == ExperimentOutcome::Reflected ||
        o == ExperimentOutcome::Escaped;
}

ExperimentOutcome classify_atmosphere_outcome(const RunSummary& r) {
    switch (r.stop_reason) {
    case StopReason::InnerBoundary:
        return ExperimentOutcome::AtmosphereReached;
    case StopReason::Mirror:
        return ExperimentOutcome::Reflected;
    case StopReason::OuterBoundary:
        return ExperimentOutcome::Escaped;
    case StopReason::TimeLimit:
        return ExperimentOutcome::TimeLimit;
    case StopReason::StepLimit:
        return ExperimentOutcome::StepLimit;
    case StopReason::NumericalError:
        return ExperimentOutcome::NumericalError;
    }
    throw std::invalid_argument("Unknown stop reason");
}

PredictionAgreement compare_prediction(std::optional<LossConePrediction> p, ExperimentOutcome o) {
    if (!is_resolved(o))
        return PredictionAgreement::Unresolved;

    if (!p || o == ExperimentOutcome::Escaped)
        return PredictionAgreement::NotApplicable;

    if (*p == LossConePrediction::Boundary)
        return PredictionAgreement::Boundary;

    const bool matches =
        (*p == LossConePrediction::AtmosphereEntry && o == ExperimentOutcome::AtmosphereReached) ||
        (*p == LossConePrediction::Reflection && o == ExperimentOutcome::Reflected);
    return matches ? PredictionAgreement::Match : PredictionAgreement::Mismatch;
}

PreparedExperiment prepare_experiment(const ExperimentParameters& p) {
    for (double value : {p.energy_keV, p.pitch_deg, p.gyro_phase_deg, p.equatorial_distance_RE,
             p.dt_s, p.max_time_s, p.atmosphere_height_m, p.outer_distance_RE,
             p.dipole_equatorial_T, p.event_tolerance_s, p.mirror_hysteresis})
        if (!std::isfinite(value))
            throw std::invalid_argument("Experiment parameters must be finite");

    if (p.energy_keV < 0 || p.pitch_deg < 0 || p.pitch_deg > 180 || p.dt_s <= 0 ||
        p.max_time_s < 0 || p.atmosphere_height_m <= 0 || p.sample_every == 0 ||
        p.event_tolerance_s <= 0 || p.mirror_hysteresis <= 0 || p.mirror_hysteresis >= 1)
        throw std::invalid_argument(
            "Invalid energy, angle, step, duration, height or event settings");

    const double atmosphere_radius = constants::earth_radius_m + p.atmosphere_height_m;
    const double initial_radius = p.equatorial_distance_RE * constants::earth_radius_m;
    const double outer_radius = p.outer_distance_RE * constants::earth_radius_m;
    if (!std::isfinite(atmosphere_radius) || !std::isfinite(initial_radius) ||
        !std::isfinite(outer_radius) ||
        !(atmosphere_radius < initial_radius && initial_radius < outer_radius) ||
        !std::isfinite(static_cast<double>(p.max_steps) * p.dt_s))
        throw std::invalid_argument(
            "Require R_atmosphere < R_initial < R_outer and a finite run duration");

    const DipoleField field(p.dipole_equatorial_T, p.dipole_axis);
    const Vec3 a = field.axis();
    const Vec3 reference = std::abs(a.x) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
    const Vec3 radial = normalized(reference - a * dot(reference, a));
    const Vec3 b = -a;

    const double pitch = p.pitch_deg * std::numbers::pi / 180;
    const double phase = std::remainder(p.gyro_phase_deg, 360.0) * std::numbers::pi / 180;
    const Vec3 direction = b * std::cos(pitch) +
        (radial * std::cos(phase) + cross(b, radial) * std::sin(phase)) * std::sin(pitch);

    const ParticleState initial{radial * initial_radius,
        momentum_from_energy(constants::eV_to_J(p.energy_keV * 1000), direction), 0};

    RunConfig config{p.dt_s, p.max_steps, p.sample_every};
    config.inner_radius_m = atmosphere_radius;
    config.outer_radius_m = outer_radius;
    config.max_duration_s = p.max_time_s;
    config.detect_mirrors = true;
    config.stop_at_first_mirror = p.stop_at_first_mirror;
    config.event_tolerance_s = p.event_tolerance_s;
    config.mirror_hysteresis = p.mirror_hysteresis;

    auto cone = equatorial_loss_cone(field, initial.position_m, atmosphere_radius);
    std::optional<LossConePrediction> prediction;
    if (p.energy_keV > 0)
        prediction = predict_equatorial_outcome(cone, pitch);

    return {field, initial, config, cone, prediction};
}

ExperimentResult run_atmosphere_experiment(
    const PreparedExperiment& e, const TrajectoryObserver& observer) {
    auto run = run_trajectory(e.initial, e.field, e.config, observer);
    const auto outcome = classify_atmosphere_outcome(run);
    const auto agreement = compare_prediction(e.prediction, outcome);

    return {std::move(run), outcome, agreement};
}
} // namespace aurora
