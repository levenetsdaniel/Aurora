#ifndef AURORA_RELATIVITY_HPP
#define AURORA_RELATIVITY_HPP

#include "../math/Vector3.hpp"
#include "Constants.hpp"

#include <cmath>
#include <stdexcept>

namespace Relativity {

inline double gammaFromMomentum(const Vector3& p, double mass_kg) {
    if (!std::isfinite(mass_kg) || mass_kg <= 0.0 || !p.isFinite())
        throw std::invalid_argument("Mass must be positive and momentum finite");

    const double mc = mass_kg * C::speedOfLight_m_s;
    return std::sqrt(1.0 + p.norm2() / (mc * mc));
}

inline Vector3 velocityFromMomentum(const Vector3& p, double mass_kg) {
    return p / (gammaFromMomentum(p, mass_kg) * mass_kg);
}

inline double kineticEnergyFromMomentum_J(const Vector3& p, double mass_kg) {
    const double gamma = gammaFromMomentum(p, mass_kg);
    return p.norm2() / (mass_kg * (gamma + 1.0));
}

inline double momentumMagnitudeFromKineticEnergy_J(double energy_J, double mass_kg) {
    if (!std::isfinite(energy_J) || energy_J < 0.0 ||
        !std::isfinite(mass_kg) || mass_kg <= 0.0)
        throw std::invalid_argument("Energy must be nonnegative and mass positive");

    const double restEnergy_J = mass_kg * C::speedOfLight_m_s * C::speedOfLight_m_s;
    return std::sqrt(energy_J * (energy_J + 2.0 * restEnergy_J)) /
        C::speedOfLight_m_s;
}

}

#endif // AURORA_RELATIVITY_HPP
