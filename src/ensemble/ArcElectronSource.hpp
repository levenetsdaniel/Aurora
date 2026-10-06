#ifndef AURORA_ARCELECTRONSOURCE_HPP
#define AURORA_ARCELECTRONSOURCE_HPP

#include "ElectronSource.hpp"

#include <cmath>
#include <cstddef>
#include <numbers>
#include <stdexcept>

// Illustrative polar-arc source. Directions are uniform in solid angle over
// the configured downward pitch interval, not a measured precipitation flux.
class ArcElectronSource final : public ElectronSource {
private:
    const MagneticField& field;
    bool auroralArc;
    double minPitch_rad;
    double maxPitch_rad;
    double physicalRate_s1;
    std::size_t markerCount;

public:
    ArcElectronSource(const MagneticField& field, bool auroralArc,
                      double minPitch_rad, double maxPitch_rad,
                      double physicalRate_s1, std::size_t markerCount)
        : field(field), auroralArc(auroralArc), minPitch_rad(minPitch_rad),
          maxPitch_rad(maxPitch_rad), physicalRate_s1(physicalRate_s1),
          markerCount(markerCount) {
        if (!std::isfinite(minPitch_rad) || !std::isfinite(maxPitch_rad) ||
            minPitch_rad < 0.0 || minPitch_rad >= maxPitch_rad ||
            maxPitch_rad > std::numbers::pi / 2.0 ||
            !std::isfinite(physicalRate_s1) || physicalRate_s1 <= 0.0 ||
            markerCount == 0 || !std::isfinite(physicalRate_s1 / markerCount) ||
            physicalRate_s1 / markerCount <= 0.0) {
            throw std::invalid_argument("ArcElectronSource: invalid configuration");
        }
    }

    InitialCondition sample(Random& random, std::uint64_t id) const override {
        const double longitude = auroralArc ?
            -0.7 + 1.4 * random.uniform01() :
            2.0 * std::numbers::pi * random.uniform01();
        const double colatitude = auroralArc ?
            0.27 + 0.05 * random.uniform01() :
            0.20 + 0.17 * random.uniform01();
        const double radius_m = C::earthRadius_m + 450000.0;
        const Vector3 position_m(
            radius_m * std::sin(colatitude) * std::cos(longitude),
            radius_m * std::sin(colatitude) * std::sin(longitude),
            radius_m * std::cos(colatitude));
        const Vector3 magnetic = field.evaluate(position_m);
        const double magnitude = magnetic.norm();
        if (!magnetic.isFinite() || !std::isfinite(magnitude) || magnitude <= 0.0) {
            throw std::runtime_error("ArcElectronSource: invalid magnetic field");
        }
        const Vector3 along = magnetic * (-1.0 / magnitude);
        const Vector3 reference = std::abs(along.z) < 0.9 ?
            Vector3(0.0, 0.0, 1.0) : Vector3(1.0, 0.0, 0.0);
        const Vector3 across = (reference - along * reference.dot(along)).normalize();
        const Vector3 around = along.cross(across);
        const double minCosine = std::cos(maxPitch_rad);
        const double maxCosine = std::cos(minPitch_rad);
        const double cosine = minCosine + (maxCosine - minCosine) * random.uniform01();
        const double pitch = std::acos(cosine);
        const double phase = 2.0 * std::numbers::pi * random.uniform01();
        const Vector3 direction = along * cosine +
            (across * std::cos(phase) + around * std::sin(phase)) * std::sin(pitch);
        const double energy_eV = 5000.0 + 20000.0 * random.uniform01();
        const double momentum = Relativity::momentumMagnitudeFromKineticEnergy_J(
            energy_eV * C::elementaryCharge_C, C::electronMass_kg);
        return {id, {position_m, direction * momentum, 0.0}, energy_eV,
                pitch, phase, physicalRate_s1 / static_cast<double>(markerCount)};
    }

    double totalPhysicalRate() const override { return physicalRate_s1; }
};

#endif // AURORA_ARCELECTRONSOURCE_HPP
