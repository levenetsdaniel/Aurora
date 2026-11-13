#pragma once

#include "aurora/fields/MagneticField.hpp"
#include "aurora/simulation/ParticleState.hpp"
#include <optional>

namespace aurora {
// Physical properties evaluated at one synchronized particle state.
// Integration, event handling and accumulation of run statistics live in simulation/.
struct ParticleProperties {
    double field_T{};
    double parallel_momentum_kg_m_s{};
    double perpendicular_momentum_kg_m_s{};
    double parallel_fraction{}; // p_parallel/|p|; zero for rest
    std::optional<double> pitch_angle_rad; // Undefined for a stationary electron.
    double magnetic_moment_J_T{}; // p_perp^2/(2*m*B), relativistic first invariant.
};
ParticleProperties particle_properties(const ParticleState& state, const MagneticField& field);
// Adiabatic estimate for analysis; the runner detects reflection from the trajectory.
double predicted_mirror_field_T(const ParticleState& initial, const MagneticField& field);
} // namespace aurora
