#pragma once

#include "aurora/math/Vec3.hpp"

namespace aurora {
// Local gyromotion scales for an electron, in SI. These do not advance a trajectory.
double gyrofrequency_rad_s(Vec3 momentum_kg_m_s, Vec3 field_T); // Non-negative; zero for B=0.
double gyroperiod_s(Vec3 momentum_kg_m_s, Vec3 field_T);         // Requires B!=0.
double larmor_radius_m(Vec3 momentum_kg_m_s, Vec3 field_T);      // Requires B!=0.
} // namespace aurora
