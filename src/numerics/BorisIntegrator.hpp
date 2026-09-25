#ifndef AURORA_BORISINTEGRATOR_HPP
#define AURORA_BORISINTEGRATOR_HPP

#include "../math/Vector3.hpp"
#include "../physics/Constants.hpp"
#include "../physics/MagneticField.hpp"

#include <cmath>
#include <stdexcept>

struct ParticleState {
    Vector3 r;
    Vector3 p;
    double time_s{};
};

class BorisIntegrator {
private:
    const MagneticField& field_;
    double charge_C_;
    double mass_kg_;

    double gammaFromMomentum(const Vector3& p) const {
        const double mc = mass_kg_ * C::speedOfLight_m_s;
        return std::sqrt(1.0 + p.norm2() / (mc * mc));
    }

    Vector3 velocityFromMomentum(const Vector3& p) const {
        return p / (gammaFromMomentum(p) * mass_kg_);
    }

public:
    BorisIntegrator(const MagneticField& field, double charge_C, double mass_kg)
        : field_(field), charge_C_(charge_C), mass_kg_(mass_kg) {
        if (!std::isfinite(charge_C_) || !std::isfinite(mass_kg_) || mass_kg_ <= 0.0)
            throw std::invalid_argument("BorisIntegrator requires a finite charge and positive finite mass");
    }

    ParticleState step(const ParticleState& state, double dt_s) const {
        if (!std::isfinite(dt_s) || dt_s <= 0.0)
            throw std::invalid_argument("BorisIntegrator step size must be positive and finite");
        if (!std::isfinite(state.time_s) || !std::isfinite(state.r.norm2()) ||
            !std::isfinite(state.p.norm2()))
            throw std::invalid_argument("ParticleState must contain finite position, momentum, and time");

        const Vector3 v0 = velocityFromMomentum(state.p);
        const Vector3 rHalf = state.r + v0 * (0.5 * dt_s);

        const Vector3 B = field_.evaluate(rHalf);
        const double gamma = gammaFromMomentum(state.p);
        const Vector3 t = B * (charge_C_ * dt_s / (2.0 * gamma * mass_kg_));
        const Vector3 s = t * (2.0 / (1.0 + t.norm2()));
        const Vector3 pPrime = state.p + state.p.cross(t);
        const Vector3 pNew = state.p + pPrime.cross(s);

        const Vector3 vNew = velocityFromMomentum(pNew);
        return {rHalf + vNew * (0.5 * dt_s), pNew, state.time_s + dt_s};
    }
};

#endif //AURORA_BORISINTEGRATOR_HPP
