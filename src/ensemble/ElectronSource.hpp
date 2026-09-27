#ifndef AURORA_ELECTRONSOURCE_HPP
#define AURORA_ELECTRONSOURCE_HPP

#include "InitialCondition.hpp"
#include "../numerics/Random.hpp"
#include "../physics/Constants.hpp"
#include "../physics/MagneticField.hpp"
#include "../physics/Relativity.hpp"

#include <cmath>
#include <cstdint>
#include <numbers>
#include <stdexcept>

class ElectronSource {
public:
    virtual ~ElectronSource() = default;
    virtual InitialCondition sample(Random& rng, std::uint64_t id) const = 0;
    virtual double totalPhysicalRate() const = 0;
};

// Uniform in energy and solid angle over one hemisphere at a fixed injection point.
// differentialFlux has units m^-2 s^-1 sr^-1 eV^-1; sourceArea is an effective m^2.
// The effective area does not add a separate cos(incidence) projection factor.
class IsotropicElectronSource final : public ElectronSource {
private:
    Vector3 position_m;
    Vector3 along;
    Vector3 transverseX;
    Vector3 transverseY;
    double time_s;
    double mass_kg;
    double minEnergy_eV;
    double maxEnergy_eV;
    double physicalRate_s1;
    double markerWeight_s1;
    std::uint64_t markersPerBatch;

public:
    IsotropicElectronSource(const MagneticField& field, Vector3 position_m,
                            double time_s, double mass_kg,
                            double minEnergy_eV, double maxEnergy_eV,
                            double differentialFlux_m2_s_sr_eV,
                            double sourceArea_m2, std::uint64_t markersPerBatch,
                            bool alongPositiveField = false)
        : position_m(position_m), time_s(time_s), mass_kg(mass_kg),
          minEnergy_eV(minEnergy_eV), maxEnergy_eV(maxEnergy_eV),
          markersPerBatch(markersPerBatch) {
        if (!position_m.isFinite() || !std::isfinite(time_s) || time_s < 0.0 ||
            !std::isfinite(mass_kg) || mass_kg <= 0.0 ||
            !std::isfinite(minEnergy_eV) || minEnergy_eV <= 0.0 ||
            !std::isfinite(maxEnergy_eV) || maxEnergy_eV <= minEnergy_eV ||
            !std::isfinite(differentialFlux_m2_s_sr_eV) ||
            differentialFlux_m2_s_sr_eV <= 0.0 ||
            !std::isfinite(sourceArea_m2) || sourceArea_m2 <= 0.0 ||
            markersPerBatch == 0) {
            throw std::invalid_argument("IsotropicElectronSource: invalid source parameters");
        }

        const Vector3 B = field.evaluate(position_m);
        const double B2 = B.norm2();
        if (!B.isFinite() || !std::isfinite(B2) || B2 <= 0.0) {
            throw std::invalid_argument("IsotropicElectronSource: field must be finite and nonzero");
        }
        const Vector3 b = B / std::sqrt(B2);
        along = alongPositiveField ? b : b * -1.0;
        const Vector3 reference = std::abs(b.x) < 0.9 ?
            Vector3(1.0, 0.0, 0.0) : Vector3(0.0, 1.0, 0.0);
        const Vector3 projected = reference - b * reference.dot(b);
        transverseX = projected / projected.norm();
        transverseY = along.cross(transverseX);

        physicalRate_s1 = differentialFlux_m2_s_sr_eV * sourceArea_m2 *
                          (maxEnergy_eV - minEnergy_eV) * (2.0 * std::numbers::pi);
        // The sampling PDF matches the uniform physical distribution, so each
        // marker in a complete batch represents the same share of the rate.
        markerWeight_s1 = physicalRate_s1 / static_cast<double>(markersPerBatch);
        if (!std::isfinite(physicalRate_s1) || physicalRate_s1 <= 0.0 ||
            !std::isfinite(markerWeight_s1) || markerWeight_s1 <= 0.0) {
            throw std::invalid_argument("IsotropicElectronSource: physical rate is not finite");
        }
    }

    InitialCondition sample(Random& rng, std::uint64_t id) const override {
        const double energy_eV = minEnergy_eV +
            (maxEnergy_eV - minEnergy_eV) * rng.uniform01();
        const double cosine = rng.uniform01();
        const double pitch_rad = std::acos(cosine);
        const double gyroPhase_rad = 2.0 * std::numbers::pi * rng.uniform01();
        const double sine = std::sqrt(1.0 - cosine * cosine);
        const Vector3 direction = along * cosine +
            (transverseX * std::cos(gyroPhase_rad) +
             transverseY * std::sin(gyroPhase_rad)) * sine;
        const double momentum = Relativity::momentumMagnitudeFromKineticEnergy_J(
            energy_eV * C::elementaryCharge_C, mass_kg);
        if (!std::isfinite(momentum) || momentum <= 0.0) {
            throw std::runtime_error("IsotropicElectronSource: sampled momentum is invalid");
        }
        return {id, {position_m, direction * momentum, time_s}, energy_eV,
                pitch_rad, gyroPhase_rad, markerWeight_s1};
    }

    double totalPhysicalRate() const override {
        return physicalRate_s1;
    }

    std::uint64_t batchSize() const {
        return markersPerBatch;
    }
};

#endif // AURORA_ELECTRONSOURCE_HPP
