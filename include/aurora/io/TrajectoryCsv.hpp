#pragma once

#include "aurora/simulation/ParticleState.hpp"
#include "aurora/fields/MagneticField.hpp"
#include <cstddef>
#include <iosfwd>

namespace aurora {
void write_trajectory_header(std::ostream& out);
void write_trajectory_row(std::ostream& out, std::size_t step, const ParticleState& state);
void write_magnetic_trajectory_header(std::ostream& out);
void write_magnetic_trajectory_row(std::ostream& out, std::size_t step,
    const ParticleState& state, const MagneticField& field);
} // namespace aurora
