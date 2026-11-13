#include "aurora/physics/LossCone.hpp"
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace {
using namespace aurora;
constexpr double R = constants::earth_radius_m;
constexpr double pi = std::numbers::pi;

void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance, const char* message) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance)
        throw std::runtime_error(message);
}

template <class F> void rejects(F action) {
    try {
        action();
    } catch (const std::exception&) {
        return;
    }
    throw std::runtime_error("Expected invalid loss cone input to be rejected");
}

void geometry(const DipoleField& field, Vec3 start, double atmosphere_radius) {
    const auto cone = equatorial_loss_cone(field, start, atmosphere_radius);
    const double initial_radius = norm(start);
    for (Vec3 footpoint : {cone.along_field_footpoint_m, cone.against_field_footpoint_m}) {
        near(
            norm(footpoint), atmosphere_radius, 1e-8, "Footpoint is not on the atmospheric sphere");
        const double sin_latitude = dot(footpoint / norm(footpoint), field.axis());
        // r_eq = r/cos^2(lambda): both footpoints must belong to the starting field line.
        near(norm(footpoint) / (1 - sin_latitude * sin_latitude), initial_radius, 1e-7,
            "Footpoint belongs to a different dipole line");
    }
    const double along_B = norm(field.at(cone.along_field_footpoint_m, 0));
    const double against_B = norm(field.at(cone.against_field_footpoint_m, 0));
    near(along_B / against_B, 1, 2e-15, "Hemisphere field magnitudes differ");
    near(cone.atmosphere_field_T, along_B, 0, "Wrong atmosphere field in result");
    near(cone.initial_field_T, norm(field.at(start, 0)), 0, "Wrong initial field in result");
    require(dot(cone.along_field_footpoint_m - start, field.at(start, 0)) > 0,
        "Along-field footpoint has wrong direction");
    require(dot(cone.against_field_footpoint_m - start, field.at(start, 0)) < 0,
        "Against-field footpoint has wrong direction");

    // Independent scalar relation: sin^2(alpha_loss) = q^3/sqrt(4-3q), q=R_atm/r_eq.
    const double q = atmosphere_radius / initial_radius;
    const double expected_ratio = q * q * q / std::sqrt(4 - 3 * q);
    near(std::pow(std::sin(cone.critical_angle_rad), 2), expected_ratio, 2e-15,
        "Critical angle differs from the scalar dipole formula");
    require(cone.critical_angle_rad > 0 && cone.critical_angle_rad < pi / 2,
        "Critical angle must be between 0 and pi/2");
}

void classification() {
    const DipoleField field;
    const auto cone = equatorial_loss_cone(field, {2 * R, 0, 0}, R + 100000);
    for (double degrees : {0., 5., 175., 180.})
        require(predict_equatorial_outcome(cone, degrees * pi / 180) ==
                LossConePrediction::AtmosphereEntry,
            "Pitch inside either cone must predict entry");
    for (double degrees : {45., 90., 135.})
        require(
            predict_equatorial_outcome(cone, degrees * pi / 180) == LossConePrediction::Reflection,
            "Pitch outside both cones must predict reflection");
    for (double angle : {cone.critical_angle_rad, pi - cone.critical_angle_rad})
        require(predict_equatorial_outcome(cone, angle) == LossConePrediction::Boundary,
            "Critical pitch must be recognized as a boundary");
    require(predict_equatorial_outcome(cone, cone.critical_angle_rad - 1e-7) ==
            LossConePrediction::AtmosphereEntry,
        "Angle below critical must be inside");
    require(predict_equatorial_outcome(cone, cone.critical_angle_rad + 1e-7) ==
            LossConePrediction::Reflection,
        "Angle above critical must be outside");
    const auto scaled = equatorial_loss_cone(DipoleField(7 * 3.12e-5), {2 * R, 0, 0}, R + 100000);
    near(scaled.critical_angle_rad, cone.critical_angle_rad, 2e-15,
        "Scaling dipole strength must not change the loss cone angle");
    std::cout << "Critical angle: " << cone.critical_angle_rad * 180 / pi
              << " deg; B_initial=" << cone.initial_field_T
              << " T; B_atmosphere=" << cone.atmosphere_field_T << " T\n";
}

void invalid_inputs() {
    const DipoleField field;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    for (double radius : {0., -1., nan, inf, 0.4 * R, 2 * R, 3 * R})
        rejects([&] {
            equatorial_loss_cone(field, {2 * R, 0, 0}, radius);
        });
    for (Vec3 position :
        {Vec3{}, Vec3{nan, 0, 0}, Vec3{inf, 0, 0}, Vec3{2 * R, 0, 0.01 * R}, Vec3{0, 0, 2 * R}})
        rejects([&] {
            equatorial_loss_cone(field, position, R + 100000);
        });
    const auto cone = equatorial_loss_cone(field, {2 * R, 0, 0}, R + 100000);
    for (double pitch : {-0.01, pi + 0.01, nan, inf})
        rejects([&] {
            predict_equatorial_outcome(cone, pitch);
        });
    for (double angle : {-0.01, pi, nan, inf}) {
        auto bad = cone;
        bad.critical_angle_rad = angle;
        rejects([&] {
            predict_equatorial_outcome(bad, 0.5);
        });
    }
}
}

int main() {
    try {
        geometry(DipoleField{}, {2 * R, 0, 0}, R + 100000);
        geometry(DipoleField{}, {0, -3 * R, 0}, R + 200000);
        const DipoleField tilted(3.12e-5, {1, 2, -3});
        const Vec3 equatorial_direction = normalized(cross(tilted.axis(), {1, 0, 0}));
        geometry(tilted, equatorial_direction * (2 * R), R + 100000);
        classification();
        invalid_inputs();
        std::cout << "PASS: loss cone geometry, both directions and invalid inputs\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
