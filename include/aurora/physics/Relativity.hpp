#pragma once

#include "aurora/math/Vec3.hpp"

namespace aurora {
// All functions describe an electron; internal units are SI.
double lorentz_factor(Vec3 momentum_kg_m_s);
Vec3 velocity_from_momentum(Vec3 momentum_kg_m_s);
double kinetic_energy_J(Vec3 momentum_kg_m_s);
Vec3 momentum_from_velocity(Vec3 velocity_m_s);
Vec3 momentum_from_energy(double kinetic_energy_J, Vec3 direction);
} // namespace aurora
