#ifndef AURORA_PRECIPITATIONCOLLECTOR_HPP
#define AURORA_PRECIPITATIONCOLLECTOR_HPP

#include "EnsembleRunner.hpp"
#include "../physics/Constants.hpp"
#include "../physics/MagneticField.hpp"
#include "../physics/Relativity.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <optional>
#include <stdexcept>

struct PrecipitatingElectron {
    std::uint64_t id{};
    ParticleState state;
    double energy_eV{};
    double incidenceCosine{};
    // Electrons per second represented by this precipitating marker.
    double weight{};
};

class PrecipitationCollector {
public:
    std::optional<PrecipitatingElectron> collect(
        const EnsembleRecord& record, const MagneticField&, double mass_kg) const {
        if (record.result.outcome != ParticleOutcome::Precipitated) {
            return std::nullopt;
        }

        const ParticleState& state = record.result.eventState;
        const double radius_m = state.r.norm();
        if (!state.r.isFinite() || !state.p.isFinite() ||
            !std::isfinite(state.time_s) ||
            !std::isfinite(radius_m) || radius_m <= 0.0 ||
            !std::isfinite(record.initial.weight) || record.initial.weight < 0.0) {
            throw std::invalid_argument("PrecipitationCollector: invalid event state or weight");
        }

        const Vector3 velocity = Relativity::velocityFromMomentum(state.p, mass_kg);
        const double speed_m_s = velocity.norm();
        const double energy_eV = Relativity::kineticEnergyFromMomentum_J(
            state.p, mass_kg) / C::elementaryCharge_C;
        if (!velocity.isFinite() || !std::isfinite(speed_m_s) || speed_m_s <= 0.0 ||
            !std::isfinite(energy_eV) || energy_eV < 0.0) {
            throw std::invalid_argument("PrecipitationCollector: invalid energy or velocity");
        }

        const double incidenceCosine = -velocity.dot(state.r / radius_m) / speed_m_s;
        if (!std::isfinite(incidenceCosine) || incidenceCosine <= 0.0 ||
            incidenceCosine > 1.0 + 1e-12) {
            throw std::invalid_argument("PrecipitationCollector: event momentum is not inward");
        }
        return PrecipitatingElectron{record.id, state, energy_eV,
            std::clamp(incidenceCosine, 0.0, 1.0), record.initial.weight};
    }
};

#endif // AURORA_PRECIPITATIONCOLLECTOR_HPP
