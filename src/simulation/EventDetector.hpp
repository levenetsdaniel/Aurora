#ifndef AURORA_EVENTDETECTOR_HPP
#define AURORA_EVENTDETECTOR_HPP


#include "Particle.hpp"
#include "../numerics/BorisIntegrator.hpp"

#include <cmath>
#include <stdexcept>

struct EventLocation {
    ParticleState state;
    double radiusResidual_m{};
    unsigned iterations{};
};

class SphereCrossingLocator {
private:
    double radius_m;
    double timeTolerance_s;
    unsigned maxIterations;

public:
    SphereCrossingLocator(double radius_m, double timeTolerance_s,
                          unsigned maxIterations = 100)
        : radius_m(radius_m), timeTolerance_s(timeTolerance_s),
          maxIterations(maxIterations) {
        if (!std::isfinite(this->radius_m) || this->radius_m <= 0.0 ||
            !std::isfinite(this->timeTolerance_s) || this->timeTolerance_s <= 0.0 ||
            this->maxIterations == 0) {
            throw std::invalid_argument("SphereCrossingLocator: invalid radius, tolerance or iteration limit");
        }
    }

    EventLocation locate(const ParticleState& before, const ParticleState& after,
                         const BorisIntegrator& integrator) const {
        if (!before.r.isFinite() || !before.p.isFinite() ||
            !after.r.isFinite() || !after.p.isFinite() ||
            !std::isfinite(before.time_s) || !std::isfinite(after.time_s) ||
            after.time_s <= before.time_s) {
            throw std::invalid_argument("SphereCrossingLocator: invalid endpoint states");
        }

        const double gBefore = before.r.norm() - radius_m;
        const double gAfter = after.r.norm() - radius_m;
        if (!std::isfinite(gBefore) || !std::isfinite(gAfter)) {
            throw std::invalid_argument("SphereCrossingLocator: non-finite boundary distance");
        }
        const bool inward = gBefore > 0.0 && gAfter <= 0.0;
        const bool outward = gBefore < 0.0 && gAfter >= 0.0;
        if (!inward && !outward) {
            throw std::invalid_argument("SphereCrossingLocator: states do not bracket a sphere crossing");
        }

        double left_s = 0.0;
        double right_s = after.time_s - before.time_s;
        ParticleState left = before;
        ParticleState right = after;
        double gLeft = gBefore;
        double gRight = gAfter;
        unsigned iterations = 0;

        while (right_s - left_s > timeTolerance_s) {
            if (iterations == maxIterations) {
                throw std::runtime_error("SphereCrossingLocator: time tolerance not reached");
            }

            const double middle_s = left_s + (right_s - left_s) * 0.5;
            if (middle_s <= left_s || middle_s >= right_s) {
                throw std::runtime_error("SphereCrossingLocator: time tolerance is below floating-point resolution");
            }

            const ParticleState middle = integrator.step(before, middle_s);
            const double gMiddle = middle.r.norm() - radius_m;
            if (!middle.r.isFinite() || !middle.p.isFinite() ||
                !std::isfinite(middle.time_s) || !std::isfinite(gMiddle)) {
                throw std::runtime_error("SphereCrossingLocator: non-finite state during localization");
            }
            iterations++;

            if (gMiddle == 0.0) {
                return {middle, 0.0, iterations};
            }
            if ((inward && gMiddle > 0.0) || (outward && gMiddle < 0.0)) {
                left_s = middle_s;
                left = middle;
                gLeft = gMiddle;
            } else {
                right_s = middle_s;
                right = middle;
                gRight = gMiddle;
            }
        }

        if (std::abs(gLeft) <= std::abs(gRight)) {
            return {left, gLeft, iterations};
        }
        return {right, gRight, iterations};
    }
};


#endif //AURORA_EVENTDETECTOR_HPP
