#include "aurora/physics/LossCone.hpp"
#include <algorithm>
#include <numbers>

namespace aurora {
EquatorialLossCone equatorial_loss_cone(
    const DipoleField& field, Vec3 initial_position, double atmosphere_radius) {
    const double equatorial_radius = norm(initial_position);
    if (!finite(initial_position) || !std::isfinite(equatorial_radius) ||
        !std::isfinite(atmosphere_radius) || atmosphere_radius <= 0 ||
        atmosphere_radius < field.minimum_radius_m() || equatorial_radius <= atmosphere_radius)
        throw std::invalid_argument(
            "Loss cone requires a finite equatorial position above the atmospheric sphere, within the field domain");

    const Vec3 radial_direction = initial_position / equatorial_radius;
    const Vec3 axis = field.axis();
    if (std::abs(dot(radial_direction, axis)) > 1e-12)
        throw std::invalid_argument(
            "Loss cone calculation supports only equatorial initial positions");

    // r = r_eq*cos^2(lambda): intersect this dipole line with r = R_atm.
    // cos(lambda_atm) = sqrt(R_atm/r_eq), sin(lambda_atm) = sqrt(1 - R_atm/r_eq).
    const double ratio = atmosphere_radius / equatorial_radius;
    const double cos_latitude = std::sqrt(ratio);
    const double sin_latitude = std::sqrt(std::max(0.0, 1.0 - ratio));
    const Vec3 equatorial_component = radial_direction * (atmosphere_radius * cos_latitude);
    const Vec3 axial_component = axis * (atmosphere_radius * sin_latitude);

    EquatorialLossCone result;
    // On the magnetic equator B points opposite to the moment axis.
    result.along_field_footpoint_m = equatorial_component - axial_component;
    result.against_field_footpoint_m = equatorial_component + axial_component;
    result.initial_field_T = norm(field.at(initial_position, 0));
    result.atmosphere_field_T = norm(field.at(result.along_field_footpoint_m, 0));
    const double field_ratio = result.initial_field_T / result.atmosphere_field_T;
    if (!std::isfinite(field_ratio) || field_ratio <= 0 || field_ratio > 1.0 + 1e-12)
        throw std::overflow_error("Loss cone field ratio is outside the supported range");

    // Conservation of |p| and mu=p_perp^2/(2mB): sin^2(alpha_loss) = B_initial/B_atm.
    result.critical_angle_rad = std::asin(std::sqrt(std::clamp(field_ratio, 0.0, 1.0)));

    return result;
}

LossConePrediction predict_equatorial_outcome(const EquatorialLossCone& cone, double pitch_angle) {
    if (!std::isfinite(pitch_angle) || pitch_angle < 0 || pitch_angle > std::numbers::pi ||
        !std::isfinite(cone.critical_angle_rad) || cone.critical_angle_rad < 0 ||
        cone.critical_angle_rad > std::numbers::pi / 2)
        throw std::invalid_argument(
            "Loss cone prediction requires pitch in [0,pi] and critical angle in [0,pi/2]");

    // The two symmetric cones are centered on +B and -B, i.e. alpha=0 and alpha=pi.
    const double angle_to_nearest_axis = std::min(pitch_angle, std::numbers::pi - pitch_angle);
    if (std::abs(angle_to_nearest_axis - cone.critical_angle_rad) <= 1e-12)
        return LossConePrediction::Boundary;
    return angle_to_nearest_axis < cone.critical_angle_rad ? LossConePrediction::AtmosphereEntry
                                                           : LossConePrediction::Reflection;
}
} // namespace aurora
