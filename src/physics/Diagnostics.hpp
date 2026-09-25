#ifndef AURORA_DIAGNOSTICS_HPP
#define AURORA_DIAGNOSTICS_HPP

#include "../math/Vector3.hpp"
#include "Constants.hpp"
#include "Relativity.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

struct ParticleDiagnostics {
    double energy_eV;
    double pitchAngle_rad;
    double magneticMoment_J_T;
    double parallelMomentum_kg_m_s;
};

inline ParticleDiagnostics evaluateDiagnostics(const Vector3& p, const Vector3& B, double mass_kg) {
    const double p2 = p.norm2();
    const double B2 = B.norm2();
    if (!p.isFinite() || !B.isFinite() ||
        !std::isfinite(p2) || p2 <= 0.0 || !std::isfinite(B2) || B2 <= 0.0)
        throw std::invalid_argument("Pitch angle and magnetic moment require nonzero finite momentum and field");

    const double Bnorm = std::sqrt(B2);
    const double pParallel = p.dot(B) / Bnorm;
    const double cosine = std::clamp(pParallel / std::sqrt(p2), -1.0, 1.0);
    const double pPerpendicular2 = p.cross(B).norm2() / B2;

    return {
        Relativity::kineticEnergyFromMomentum_J(p, mass_kg) / C::elementaryCharge_C,
        std::acos(cosine),
        pPerpendicular2 / (2.0 * mass_kg * Bnorm),
        pParallel
    };
}

#endif // AURORA_DIAGNOSTICS_HPP
