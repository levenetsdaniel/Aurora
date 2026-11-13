#include "aurora/simulation/TrajectoryRunner.hpp"
#include "aurora/simulation/EventLocator.hpp"
#include "aurora/physics/ParticleProperties.hpp"
#include "aurora/physics/Relativity.hpp"
#include <algorithm>

namespace aurora {
const char* stop_reason_name(StopReason reason) {
    switch (reason) {
    case StopReason::StepLimit:
        return "step_limit";
    case StopReason::TimeLimit:
        return "time_limit";
    case StopReason::InnerBoundary:
        return "inner_boundary";
    case StopReason::OuterBoundary:
        return "outer_boundary";
    case StopReason::Mirror:
        return "mirror";
    case StopReason::NumericalError:
        return "numerical_error";
    }
    return "unknown";
}

namespace {
void validate_config(const ParticleState& initial, const RunConfig& c) {
    validate_state(initial);
    if (c.sample_every == 0 || !(c.event_tolerance_s > 0) || !std::isfinite(c.event_tolerance_s) ||
        !(c.mirror_hysteresis > 0 && c.mirror_hysteresis < 1) || !(c.max_duration_s >= 0) ||
        (c.stop_at_first_mirror && !c.detect_mirrors))
        throw std::invalid_argument("Invalid trajectory limits or event settings");

    for (const auto radius : {c.inner_radius_m, c.outer_radius_m})
        if (radius && (!std::isfinite(*radius) || !(*radius > 0)))
            throw std::invalid_argument("Boundary radii must be finite and positive");
    if (c.inner_radius_m && c.outer_radius_m && *c.inner_radius_m >= *c.outer_radius_m)
        throw std::invalid_argument("Inner radius must be less than outer radius");

    const double r = norm(initial.position_m);
    if ((c.inner_radius_m && r <= *c.inner_radius_m) ||
        (c.outer_radius_m && r >= *c.outer_radius_m))
        throw std::invalid_argument(
            "Initial particle must be strictly inside the spherical domain");

    const double requested = initial.time_s + static_cast<double>(c.steps) * c.dt_s;
    if (!std::isfinite(requested) || (c.steps > 0 && !(requested > initial.time_s)))
        throw std::invalid_argument("Requested duration is outside the supported range");
    if (std::isfinite(c.max_duration_s) && !std::isfinite(initial.time_s + c.max_duration_s))
        throw std::invalid_argument("Time limit overflow");
}

int sign(double value) {
    return value > 0 ? 1 : (value < 0 ? -1 : 0);
}
}

RunSummary run_trajectory(const ParticleState& initial, const MagneticField& field,
    RunConfig config, const TrajectoryObserver& observer) {
    const RelativisticBorisPusher pusher(config.dt_s);
    validate_config(initial, config);
    auto state = pusher.initialize(initial, field);

    RunSummary summary;
    summary.final_state = initial;
    summary.initial_energy_J = kinetic_energy_J(initial.momentum_kg_m_s);
    summary.max_speed_m_s = norm(velocity_from_momentum(initial.momentum_kg_m_s));

    int armed_sign = 0, expected_sign = 0;

    if (config.detect_mirrors) {
        const auto d = particle_properties(initial, field);
        summary.initial_magnetic_moment_J_T = d.magnetic_moment_J_T;
        summary.max_field_T = d.field_T;
        const double u = d.parallel_fraction;
        if (std::abs(u) > config.mirror_hysteresis)
            armed_sign = sign(u);
    }

    if (observer)
        observer(0, initial);

    if (config.max_duration_s == 0) {
        summary.stop_reason = StopReason::TimeLimit;

        return summary;
    }

    const double time_limit = initial.time_s + config.max_duration_s;

    for (std::size_t i = 0; i < config.steps; ++i) {
        const ParticleState before = summary.final_state;
        ParticleState sample;
        bool terminal = false;
        bool partial = false;

        try {
            const double full_end = state.time_s() + config.dt_s;
            partial = time_limit <= full_end;
            const double h = partial ? time_limit - before.time_s : config.dt_s;
            if (!(h > 0)) {
                summary.stop_reason = StopReason::TimeLimit;
                break;
            }

            auto next = state;
            Vec3 end_position;
            if (partial) {
                const RelativisticBorisPusher last_step(h);
                next = last_step.initialize(before, field);
                last_step.advance(next, field);
            } else
                pusher.advance(next, field);
            end_position = next.position_m();

            std::optional<TrajectoryEvent> boundary;
            if (config.inner_radius_m) {
                const double hi =
                    inner_crossing_bracket(before.position_m, end_position, *config.inner_radius_m);
                if (hi > 0) {
                    const auto e = locate_radius_event(
                        before, field, h, *config.inner_radius_m, 0, hi, config.event_tolerance_s);
                    boundary =
                        TrajectoryEvent{EventType::InnerBoundary, e.state, e.bracket_width_s};
                }
            }

            if (config.outer_radius_m && norm(end_position) >= *config.outer_radius_m) {
                const auto e = locate_radius_event(
                    before, field, h, *config.outer_radius_m, 0, 1, config.event_tolerance_s);
                if (!boundary || e.state.time_s < boundary->state.time_s)
                    boundary =
                        TrajectoryEvent{EventType::OuterBoundary, e.state, e.bracket_width_s};
            }

            // Boundary detection precedes endpoint field evaluation: overshooting the field cutoff is safe.
            sample =
                boundary ? boundary->state : RelativisticBorisPusher(h).synchronize(next, field);

            if (config.detect_mirrors) {
                const double u0 = particle_properties(before, field).parallel_fraction;
                const double u1 = particle_properties(sample, field).parallel_fraction;
                // Для выбранного направления s ищем переход s*u_0 > 0 к s*u_1 <= 0: разворот.
                if (armed_sign != 0 && u0 * armed_sign > 0 && u1 * armed_sign <= 0) {
                    const double mirror_interval = sample.time_s - before.time_s;
                    const auto e = locate_mirror_event(
                        before, field, mirror_interval, config.event_tolerance_s);
                    summary.events.push_back({EventType::Mirror, e.state, e.bracket_width_s});
                    expected_sign = -armed_sign;
                    armed_sign = 0;
                    if (config.stop_at_first_mirror) {
                        sample = e.state;
                        boundary.reset();
                        summary.stop_reason = StopReason::Mirror;
                        terminal = true;
                    }
                }
                // Повторное включение: |u_1| > hysteresis и sgn(u_1) = ожидаемому направлению.
                if (armed_sign == 0 && std::abs(u1) > config.mirror_hysteresis &&
                    (expected_sign == 0 || sign(u1) == expected_sign))
                    armed_sign = sign(u1);
            }

            if (boundary) {
                summary.events.push_back(*boundary);
                summary.stop_reason = boundary->type == EventType::InnerBoundary
                    ? StopReason::InnerBoundary
                    : StopReason::OuterBoundary;
                terminal = true;
            } else if (partial && !terminal) {
                summary.stop_reason = StopReason::TimeLimit;
                terminal = true;
            }

            if (!terminal)
                state = next;

            if (config.detect_mirrors) {
                const auto d = particle_properties(sample, field);
                summary.max_field_T = std::max(summary.max_field_T, d.field_T);
                if (summary.initial_magnetic_moment_J_T > 0)
                    summary.max_relative_magnetic_moment_error = std::max(
                        summary.max_relative_magnetic_moment_error,
                        std::abs(d.magnetic_moment_J_T / summary.initial_magnetic_moment_J_T - 1));
            }

            const double energy = kinetic_energy_J(sample.momentum_kg_m_s);
            if (summary.initial_energy_J > 0)
                summary.max_relative_energy_error = std::max(summary.max_relative_energy_error,
                    std::abs(energy / summary.initial_energy_J - 1));
            summary.max_speed_m_s = std::max(
                summary.max_speed_m_s, norm(velocity_from_momentum(sample.momentum_kg_m_s)));
        } catch (const std::exception& e) {
            summary.stop_reason = StopReason::NumericalError;
            summary.error_message = e.what();
            // The last valid state is the result. Observer/I/O errors are not swallowed.
            if (observer && i % config.sample_every != 0)
                observer(i, summary.final_state);

            return summary;
        }

        summary.completed_steps = i + 1; // A localized final fraction counts as an attempted step.
        summary.final_state = sample;
        if (observer && ((i + 1) % config.sample_every == 0 || i + 1 == config.steps || terminal))
            observer(i + 1, sample);

        if (terminal)
            break;
    }

    return summary;
}
} // namespace aurora
