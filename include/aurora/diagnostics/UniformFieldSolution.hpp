#pragma once

#include "aurora/simulation/ParticleState.hpp"

namespace aurora {
// Independent continuous-time reference; NOT used by the numerical pusher.
ParticleState uniform_field_solution(const ParticleState& initial, Vec3 field_T, double time_s);
double gyrofrequency_rad_s(Vec3 momentum, Vec3 field_T); // Positive magnitude.
double gyroperiod_s(Vec3 momentum, Vec3 field_T);         // Rejects B=0.
double larmor_radius_m(Vec3 momentum, Vec3 field_T);      // Rejects B=0.
} // namespace aurora
