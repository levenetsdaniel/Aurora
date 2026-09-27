#ifndef AURORA_LOSSCONE_HPP
#define AURORA_LOSSCONE_HPP


#include "SimulationEngine.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>

struct LossConeBracket {
    double precipitatingPitch_rad{};
    double mirroringPitch_rad{};
    unsigned iterations{};
};

class LossConeScanner {
private:
    const SimulationEngine& engine;
    unsigned maxIterations;

    ParticleOutcome outcomeAtPitch(const ParticleState& source, const Vector3& along,
                                   const Vector3& perpendicular, double momentum,
                                   double pitch_rad) const {
        ParticleState trial = source;
        trial.p = (along * std::cos(pitch_rad) +
                   perpendicular * std::sin(pitch_rad)) * momentum;
        return engine.run(trial).outcome;
    }

public:
    explicit LossConeScanner(const SimulationEngine& engine, unsigned maxIterations = 100)
        : engine(engine), maxIterations(maxIterations) {
        if (maxIterations == 0) {
            throw std::invalid_argument("LossConeScanner: maxIterations must be positive");
        }
    }

    LossConeBracket refine(const ParticleState& equatorialTemplate,
                           double lowPitch_rad, double highPitch_rad,
                           double tolerance_rad) const {
        if (!std::isfinite(lowPitch_rad) || !std::isfinite(highPitch_rad) ||
            !std::isfinite(tolerance_rad) || lowPitch_rad < 0.0 ||
            lowPitch_rad >= highPitch_rad ||
            highPitch_rad >= std::numbers::pi / 2.0 || tolerance_rad <= 0.0) {
            throw std::invalid_argument("LossConeScanner: invalid pitch bracket or tolerance");
        }
        if (!equatorialTemplate.r.isFinite() || !equatorialTemplate.p.isFinite() ||
            !std::isfinite(equatorialTemplate.time_s)) {
            throw std::invalid_argument("LossConeScanner: invalid template state");
        }

        const Vector3 B = engine.magneticField().evaluate(equatorialTemplate.r);
        const double B2 = B.norm2();
        const double momentum2 = equatorialTemplate.p.norm2();
        if (!B.isFinite() || !std::isfinite(B2) || B2 <= 0.0 ||
            !std::isfinite(momentum2) || momentum2 <= 0.0) {
            throw std::invalid_argument("LossConeScanner: field and momentum must be finite and nonzero");
        }

        const Vector3 b = B / std::sqrt(B2);
        const double parallel = equatorialTemplate.p.dot(b);
        const Vector3 transverse = equatorialTemplate.p - b * parallel;
        const double transverseNorm = transverse.norm();
        if (!std::isfinite(parallel) || parallel == 0.0 ||
            !std::isfinite(transverseNorm) || transverseNorm <= 0.0) {
            throw std::invalid_argument("LossConeScanner: template must define a hemisphere and gyro phase");
        }

        const Vector3 along = parallel > 0.0 ? b : b * -1.0;
        const Vector3 perpendicular = transverse / transverseNorm;
        const double momentum = std::sqrt(momentum2);

        if (outcomeAtPitch(equatorialTemplate, along, perpendicular, momentum,
                           lowPitch_rad) != ParticleOutcome::Precipitated) {
            throw std::invalid_argument("LossConeScanner: lower pitch must precipitate");
        }
        if (outcomeAtPitch(equatorialTemplate, along, perpendicular, momentum,
                           highPitch_rad) != ParticleOutcome::Mirrored) {
            throw std::invalid_argument("LossConeScanner: upper pitch must mirror");
        }

        unsigned iterations = 0;
        while (highPitch_rad - lowPitch_rad > tolerance_rad) {
            if (iterations == maxIterations) {
                throw std::runtime_error("LossConeScanner: pitch tolerance not reached");
            }

            const double middle = lowPitch_rad + (highPitch_rad - lowPitch_rad) * 0.5;
            if (middle <= lowPitch_rad || middle >= highPitch_rad) {
                throw std::runtime_error("LossConeScanner: pitch tolerance is below floating-point resolution");
            }

            const ParticleOutcome outcome = outcomeAtPitch(
                equatorialTemplate, along, perpendicular, momentum, middle);
            if (outcome == ParticleOutcome::Precipitated) {
                lowPitch_rad = middle;
            } else if (outcome == ParticleOutcome::Mirrored) {
                highPitch_rad = middle;
            } else {
                throw std::runtime_error("LossConeScanner: midpoint escaped or remained unresolved");
            }
            iterations++;
        }

        return {lowPitch_rad, highPitch_rad, iterations};
    }
};


#endif // AURORA_LOSSCONE_HPP
