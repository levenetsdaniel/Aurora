#ifndef AURORA_BORISINTEGRATOR_HPP
#define AURORA_BORISINTEGRATOR_HPP

#include "../math/Vector3.hpp"
#include "../physics/MagneticField.hpp"
#include "../physics/Relativity.hpp"
#include "../simulation/Particle.hpp"

#include <cmath>
#include <stdexcept>


class BorisIntegrator {
private:
    const MagneticField& field;
    double charge_C;
    double mass_kg;

public:
    BorisIntegrator(const MagneticField& field, double charge_C, double mass_kg)
        : field(field), charge_C(charge_C), mass_kg(mass_kg) {
        if (!std::isfinite(this->charge_C) || !std::isfinite(this->mass_kg) || this->mass_kg <= 0.0)
            throw std::invalid_argument("BorisIntegrator requires a finite charge and positive finite mass");
    }

    ParticleState step(const ParticleState& state, double dt_s) const {
        if (!std::isfinite(dt_s) || dt_s <= 0.0)
            throw std::invalid_argument("BorisIntegrator step size must be positive and finite");
        if (!std::isfinite(state.time_s) || !state.r.isFinite() || !state.p.isFinite())
            throw std::invalid_argument("ParticleState must contain finite position, momentum, and time");

        const Vector3 v0 = Relativity::velocityFromMomentum(state.p, mass_kg);
        const Vector3 rHalf = state.r + v0 * (0.5 * dt_s);

        const Vector3 B = field.evaluate(rHalf);
        const double gamma = Relativity::gammaFromMomentum(state.p, mass_kg);
        const Vector3 t = B * (charge_C * dt_s / (2.0 * gamma * mass_kg));
        const Vector3 s = t * (2.0 / (1.0 + t.norm2()));
        const Vector3 pPrime = state.p + state.p.cross(t);
        const Vector3 pNew = state.p + pPrime.cross(s);

        const Vector3 vNew = Relativity::velocityFromMomentum(pNew, mass_kg);
        return {rHalf + vNew * (0.5 * dt_s), pNew, state.time_s + dt_s};
    }
};

#endif //AURORA_BORISINTEGRATOR_HPP
