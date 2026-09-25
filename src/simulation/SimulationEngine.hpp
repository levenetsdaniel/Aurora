#ifndef AURORA_SIMULATIONENGINE_HPP
#define AURORA_SIMULATIONENGINE_HPP


#include "../physics/MagneticField.hpp"
#include "../physics/Constants.hpp"
#include "../numerics/BorisIntegrator.hpp"
#include "SimulationConfig.hpp"
#include "SimulationResult.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

class SimulationEngine {
private:
    const MagneticField& field;
    double charge_C;
    double mass_kg;
    SimulationConfig config;

    void validateConfig() const {
        if (!std::isfinite(config.dt_s) || config.dt_s <= 0.0) {
            throw std::invalid_argument("SimulationConfig: dt_s must be finite and > 0");
        }

        if (!std::isfinite(config.atmosphereHeight_m) ||
            config.atmosphereHeight_m < 0.0) {
            throw std::invalid_argument("SimulationConfig: atmosphereHeight_m must be finite and >= 0");
            }

        if (!std::isfinite(config.escapeRadius_m) ||
            config.escapeRadius_m <= 0.0) {
            throw std::invalid_argument("SimulationConfig: escapeRadius_m must be finite and > 0");
            }

        const double atmosphereRadius_m =
            C::earthRadius_m + config.atmosphereHeight_m;

        if (config.escapeRadius_m <= atmosphereRadius_m) {
            throw std::invalid_argument("SimulationConfig: escapeRadius_m must be greater than atmosphere radius");
        }

        if (!std::isfinite(config.maxTime_s) ||
        config.maxTime_s <= 0.0) {
            throw std::invalid_argument("SimulationConfig: maxTime_s must be finite and > 0");
        }

        if (config.maxSteps == 0) {
            throw std::invalid_argument("SimulationConfig: maxSteps must be > 0");
        }
    }

    void validateState(const ParticleState& state) const {
        if (!state.r.isFinite()) {
            throw std::invalid_argument("ParticleState: position must contain only finite values");
        }

        if (!state.p.isFinite()) {
            throw std::invalid_argument("ParticleState: momentum must contain only finite values");
        }

        if (!std::isfinite(state.time_s)) {
            throw std::invalid_argument("ParticleState: time_s must be finite");
        }

        if (state.time_s < 0.0) {
            throw std::invalid_argument("ParticleState: time_s must be >= 0");
        }

        const double radius_m = state.r.norm();
        const double atmosphereRadius_m = C::earthRadius_m + config.atmosphereHeight_m;
        if (!std::isfinite(radius_m) || radius_m <= atmosphereRadius_m ||
            radius_m >= config.escapeRadius_m) {
            throw std::invalid_argument("ParticleState: initial position must be between the two boundaries");
        }

        if (state.time_s >= config.maxTime_s) {
            throw std::invalid_argument("ParticleState: initial time must be below maxTime_s");
        }
    }

public:
    SimulationEngine(const MagneticField& field, double chargeC, double mass_kg, SimulationConfig config)
        : field(field), charge_C(chargeC), mass_kg(mass_kg), config(config) {
        validateConfig();
    }

    SimulationResult run(const ParticleState& initial) const {
        validateState(initial);

        const double atmosphereRadius_m = config.atmosphereHeight_m + C::earthRadius_m;

        BorisIntegrator boris(field, charge_C, mass_kg);
        ParticleState current = initial;
        std::uint64_t steps = 0;

        SimulationResult result{};

        while (steps < config.maxSteps && current.time_s < config.maxTime_s) {
            const double remaining_s = config.maxTime_s - current.time_s;
            if (remaining_s <= 0.0) break;

            const double step_s = std::min(config.dt_s, remaining_s);
            const ParticleState next = boris.step(current, step_s);
            steps++;

            const double rCurrent_m = current.r.norm();
            const double rNext_m = next.r.norm();

            if (rCurrent_m > atmosphereRadius_m && rNext_m <= atmosphereRadius_m) {
                result.outcome = ParticleOutcome::Precipitated;
                result.eventState = next;
                result.steps = steps;
                return result;
            }

            if (rCurrent_m < config.escapeRadius_m && rNext_m >= config.escapeRadius_m) {
                result.outcome = ParticleOutcome::Escaped;
                result.eventState = next;
                result.steps = steps;
                return result;
            }

            current = next;
        }

        result.outcome = ParticleOutcome::Unresolved;
        result.eventState = current;
        result.steps = steps;
        return result;
    }
};


#endif //AURORA_SIMULATIONENGINE_HPP
