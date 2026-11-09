#include "aurora/simulation/TrajectoryRunner.hpp"
#include "aurora/physics/Relativity.hpp"
#include <algorithm>

namespace aurora {
RunSummary run_trajectory(const ParticleState& initial, const MagneticField& field,
                          RunConfig config, const TrajectoryObserver& observer) {
    if (config.sample_every == 0) throw std::invalid_argument("sample_every must be positive");
    RelativisticBorisPusher pusher(config.dt_s);
    auto state = pusher.initialize(initial, field);
    const double end_time = initial.time_s + static_cast<double>(config.steps)*config.dt_s;
    if (!std::isfinite(end_time) || (config.steps > 0 && !(end_time > initial.time_s)))
        throw std::invalid_argument("Requested duration is outside the supported range");
    RunSummary summary{initial, 0, kinetic_energy_J(initial.momentum_kg_m_s), 0.0,
                       norm(velocity_from_momentum(initial.momentum_kg_m_s))};
    if (observer) observer(0, initial);
    // Loop avoids wraparound even when steps == size_t::max().
    for (std::size_t i = 0; i < config.steps; ++i) {
        pusher.advance(state, field);
        const auto sample = pusher.synchronize(state, field);
        const double energy = kinetic_energy_J(sample.momentum_kg_m_s);
        if (summary.initial_energy_J > 0.0)
            summary.max_relative_energy_error = std::max(summary.max_relative_energy_error,
                std::abs(energy/summary.initial_energy_J - 1.0));
        summary.max_speed_m_s = std::max(summary.max_speed_m_s,
            norm(velocity_from_momentum(sample.momentum_kg_m_s)));
        summary.completed_steps = i+1;
        summary.final_state = sample;
        if (observer && ((i+1)%config.sample_every == 0 || i+1 == config.steps))
            observer(i+1, sample);
    }
    return summary;
}
} // namespace aurora
