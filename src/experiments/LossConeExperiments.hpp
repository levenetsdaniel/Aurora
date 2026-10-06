#ifndef AURORA_LOSSCONEEXPERIMENTS_HPP
#define AURORA_LOSSCONEEXPERIMENTS_HPP

#include "TrajectoryExperiments.hpp"
#include "../simulation/LossCone.hpp"

#include <cmath>
#include <cstdint>
#include <numbers>
#include <stdexcept>

namespace LossConeExperiments {

struct Result {
    double sourceColatitude_deg;
    double L;
    double atmosphereHeight_m;
    int stepsPerGyroperiod;
    double dt_s;
    double precipitatingPitch_deg;
    double mirroringPitch_deg;
    double numericalPitch_deg;
    double theoreticalPitch_deg;
    double initialField_T;
    double boundaryField_T;
    unsigned iterations;
};

inline Result run(double colatitude_rad, double atmosphereHeight_m,
                  int stepsPerGyroperiod, double tolerance_deg = 0.05) {
    if (!std::isfinite(colatitude_rad) || colatitude_rad <= 0.0 ||
        colatitude_rad >= std::numbers::pi / 2.0 ||
        !std::isfinite(atmosphereHeight_m) || atmosphereHeight_m < 0.0 ||
        atmosphereHeight_m >= 450000.0 || stepsPerGyroperiod < 4 ||
        !std::isfinite(tolerance_deg) || tolerance_deg <= 0.0 ||
        tolerance_deg >= 5.0) {
        throw std::invalid_argument("LossConeExperiments: invalid configuration");
    }
    const DipoleField field;
    const double sourceRadius_m = C::earthRadius_m + 450000.0;
    const double boundaryRadius_m = C::earthRadius_m + atmosphereHeight_m;
    const double sine = std::sin(colatitude_rad);
    const double L = sourceRadius_m / (C::earthRadius_m * sine * sine);
    const double boundarySine2 = boundaryRadius_m / (L * C::earthRadius_m);
    if (boundarySine2 <= 0.0 || boundarySine2 >= 1.0)
        throw std::invalid_argument("LossConeExperiments: field line misses boundary");
    const double boundaryColatitude = std::asin(std::sqrt(boundarySine2));
    const Vector3 sourcePosition(sourceRadius_m * sine, 0.0,
                                 sourceRadius_m * std::cos(colatitude_rad));
    const Vector3 boundaryPosition(boundaryRadius_m * std::sin(boundaryColatitude),
                                   0.0, boundaryRadius_m * std::cos(boundaryColatitude));
    const Vector3 sourceField = field.evaluate(sourcePosition);
    const double initialField_T = sourceField.norm();
    const double boundaryField_T = field.evaluate(boundaryPosition).norm();
    const double theoreticalPitch_rad = std::asin(std::sqrt(
        initialField_T / boundaryField_T));
    const Vector3 towardEarth = sourceField * (-1.0 / initialField_T);
    const Vector3 transverse(0.0, 1.0, 0.0);
    constexpr double templatePitch_rad = std::numbers::pi / 4.0;
    const ParticleState initial{sourcePosition,
        (towardEarth * std::cos(templatePitch_rad) +
         transverse * std::sin(templatePitch_rad)) *
        TrajectoryExperiments::momentumMagnitude(), 0.0};
    const double dt_s = TrajectoryExperiments::gyroPeriod(initialField_T) /
        static_cast<double>(stepsPerGyroperiod);
    const SimulationEngine engine(field, -C::elementaryCharge_C, C::electronMass_kg,
        {dt_s, atmosphereHeight_m, C::earthRadius_m + 800000.0,
         0.05, static_cast<std::uint64_t>(std::ceil(0.05 / dt_s)) + 1});
    const LossConeBracket bracket = LossConeScanner(engine).refine(
        initial, 50.0 * std::numbers::pi / 180.0,
        85.0 * std::numbers::pi / 180.0,
        tolerance_deg * std::numbers::pi / 180.0);
    const double degrees = 180.0 / std::numbers::pi;
    return {colatitude_rad * degrees, L, atmosphereHeight_m,
            stepsPerGyroperiod, dt_s,
            bracket.precipitatingPitch_rad * degrees,
            bracket.mirroringPitch_rad * degrees,
            (bracket.precipitatingPitch_rad + bracket.mirroringPitch_rad) *
                0.5 * degrees,
            theoreticalPitch_rad * degrees, initialField_T,
            boundaryField_T, bracket.iterations};
}

} // namespace LossConeExperiments

#endif // AURORA_LOSSCONEEXPERIMENTS_HPP
