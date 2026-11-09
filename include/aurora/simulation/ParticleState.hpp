#pragma once

#include "aurora/math/Vec3.hpp"

namespace aurora {
struct ParticleState {
    Vec3 position_m;
    Vec3 momentum_kg_m_s;
    double time_s{};
}; // Position and momentum both refer to time_s.

inline void validate_state(const ParticleState& state) {
    if (!finite(state.position_m) || !finite(state.momentum_kg_m_s) || !std::isfinite(state.time_s))
        throw std::invalid_argument("Particle state must be finite");
}
} // namespace aurora
