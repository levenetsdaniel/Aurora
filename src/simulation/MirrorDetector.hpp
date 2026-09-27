#ifndef AURORA_MIRRORDETECTOR_HPP
#define AURORA_MIRRORDETECTOR_HPP


#include "Particle.hpp"
#include "../numerics/BorisIntegrator.hpp"
#include "../physics/MagneticField.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

struct MirrorEvent {
    ParticleState state;
    double latitude_rad{};
};

class MirrorDetector {
private:
    double parallelMomentumTolerance;
    double timeTolerance_s;
    unsigned maxIterations;

    static double parallelMomentum(const ParticleState& state, const MagneticField& field) {
        if (!state.r.isFinite() || !state.p.isFinite() || !std::isfinite(state.time_s)) {
            throw std::invalid_argument("MirrorDetector: invalid particle state");
        }

        const Vector3 B = field.evaluate(state.r);
        const double B2 = B.norm2();
        if (!B.isFinite() || !std::isfinite(B2) || B2 <= 0.0) {
            throw std::invalid_argument("MirrorDetector: magnetic field must be finite and nonzero");
        }

        const double value = state.p.dot(B) / std::sqrt(B2);
        if (!std::isfinite(value)) {
            throw std::invalid_argument("MirrorDetector: non-finite parallel momentum");
        }
        return value;
    }

    static MirrorEvent makeEvent(const ParticleState& state) {
        const double radius = state.r.norm();
        if (!state.r.isFinite() || !std::isfinite(radius) || radius <= 0.0) {
            throw std::runtime_error("MirrorDetector: latitude is undefined at the origin");
        }
        return {state, std::asin(std::clamp(state.r.z / radius, -1.0, 1.0))};
    }

public:
    explicit MirrorDetector(double parallelMomentumTolerance = 0.0,
                            double timeTolerance_s = 1e-12,
                            unsigned maxIterations = 100)
        : parallelMomentumTolerance(parallelMomentumTolerance),
          timeTolerance_s(timeTolerance_s), maxIterations(maxIterations) {
        if (!std::isfinite(this->parallelMomentumTolerance) ||
            this->parallelMomentumTolerance < 0.0 ||
            !std::isfinite(this->timeTolerance_s) || this->timeTolerance_s <= 0.0 ||
            this->maxIterations == 0) {
            throw std::invalid_argument("MirrorDetector: invalid tolerance or iteration limit");
        }
    }

    bool bracketed(const ParticleState& before, const ParticleState& after,
                   const MagneticField& field) const {
        const double beforeParallel = parallelMomentum(before, field);
        const double afterParallel = parallelMomentum(after, field);
        return (beforeParallel > parallelMomentumTolerance &&
                afterParallel < -parallelMomentumTolerance) ||
               (beforeParallel < -parallelMomentumTolerance &&
                afterParallel > parallelMomentumTolerance);
    }

    MirrorEvent locate(const ParticleState& before, const ParticleState& after,
                       const MagneticField& field, const BorisIntegrator& integrator) const {
        if (!std::isfinite(before.time_s) || !std::isfinite(after.time_s) ||
            after.time_s <= before.time_s || !bracketed(before, after, field)) {
            throw std::invalid_argument("MirrorDetector: states do not bracket a mirror event");
        }

        const double initialParallel = parallelMomentum(before, field);
        double left_s = 0.0;
        double right_s = after.time_s - before.time_s;
        unsigned iterations = 0;

        while (right_s - left_s > timeTolerance_s) {
            if (iterations == maxIterations) {
                throw std::runtime_error("MirrorDetector: time tolerance not reached");
            }

            const double middle_s = left_s + (right_s - left_s) * 0.5;
            if (middle_s <= left_s || middle_s >= right_s) {
                throw std::runtime_error("MirrorDetector: time tolerance is below floating-point resolution");
            }

            const ParticleState middle = integrator.step(before, middle_s);
            const double middleParallel = parallelMomentum(middle, field);
            if (middleParallel == 0.0) return makeEvent(middle);

            if ((initialParallel > 0.0 && middleParallel > 0.0) ||
                (initialParallel < 0.0 && middleParallel < 0.0)) {
                left_s = middle_s;
            } else {
                right_s = middle_s;
            }
            iterations++;
        }

        const double eventTime_s = left_s + (right_s - left_s) * 0.5;
        if (eventTime_s <= 0.0) {
            throw std::runtime_error("MirrorDetector: mirror time cannot be represented");
        }
        return makeEvent(integrator.step(before, eventTime_s));
    }

};


#endif //AURORA_MIRRORDETECTOR_HPP
