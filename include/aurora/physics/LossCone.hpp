#pragma once

#include "aurora/fields/DipoleField.hpp"

namespace aurora {
struct EquatorialLossCone {
    Vec3 along_field_footpoint_m;
    Vec3 against_field_footpoint_m;
    double initial_field_T{};
    double atmosphere_field_T{}; // Same magnitude at both ends of a centered dipole line.
    double critical_angle_rad{}; // Half-angle of each of the two loss cones.
};

// Analytical field-line geometry for an equatorial start in a centered dipole.
// A tilted moment axis is supported; off-equatorial starting positions are rejected.
EquatorialLossCone equatorial_loss_cone(
    const DipoleField& field, Vec3 equatorial_position_m, double atmosphere_radius_m);

enum class LossConePrediction { AtmosphereEntry, Reflection, Boundary };
// Adiabatic prediction for a moving particle; the runner determines the actual outcome.
// The boundary is recognized to within 1e-12 rad to allow for floating-point roundoff.
LossConePrediction predict_equatorial_outcome(
    const EquatorialLossCone& cone, double pitch_angle_rad);
} // namespace aurora
