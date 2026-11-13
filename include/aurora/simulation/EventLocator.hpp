#pragma once

#include "aurora/fields/MagneticField.hpp"
#include "aurora/simulation/ParticleState.hpp"

namespace aurora {
struct LocatedEvent {
    ParticleState state;
    double bracket_width_s{};
};

// A fresh fractional Boris step from a synchronized state; only used at events/time limits.
ParticleState fractional_step(
    const ParticleState& start, const MagneticField& field, double interval_s);
// Fractions bracket a sign change on one full step of length dt_s.
LocatedEvent locate_radius_event(const ParticleState& start, const MagneticField& field,
    double dt_s, double radius_m, double lower_fraction, double upper_fraction, double tolerance_s);
LocatedEvent locate_mirror_event(
    const ParticleState& start, const MagneticField& field, double dt_s, double tolerance_s);
// Finds the first inward sphere crossing, including an enter/exit within one straight drift.
// Returns a bracket upper fraction, or 0 if no crossing candidate exists.
double inner_crossing_bracket(Vec3 start, Vec3 end, double radius_m);
} // namespace aurora
