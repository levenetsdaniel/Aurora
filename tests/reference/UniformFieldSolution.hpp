#pragma once

#include "aurora/simulation/ParticleState.hpp"

namespace aurora::test_reference {
// Continuous-time analytical trajectory for comparison with the numerical integrator.
// Compiled only when BUILD_TESTING=ON; never used to advance the simulation.
ParticleState uniform_field_solution(const ParticleState& initial, Vec3 field_T, double time_s);
} // namespace aurora::test_reference
