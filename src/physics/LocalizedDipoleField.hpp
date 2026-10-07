#ifndef AURORA_LOCALIZEDDIPOLEFIELD_HPP
#define AURORA_LOCALIZEDDIPOLEFIELD_HPP

#include "Constants.hpp"
#include "DipoleField.hpp"

#include <cmath>
#include <numbers>
#include <stdexcept>

// A divergence-free vacuum field outside a dipole buried below the surface.
// surfacePeak_T is the signed axial anomaly at the surface directly above it.
class LocalizedDipoleField final : public MagneticField {
private:
    DipoleField background;
    Vector3 source;
    Vector3 axis;
    double coefficient_T_m3;

public:
    LocalizedDipoleField(double colatitude_rad, double depth_m,
                         double surfacePeak_T)
        : source((C::earthRadius_m - depth_m) * std::sin(colatitude_rad), 0.0,
                 (C::earthRadius_m - depth_m) * std::cos(colatitude_rad)),
          axis(std::sin(colatitude_rad), 0.0, std::cos(colatitude_rad)),
          coefficient_T_m3(surfacePeak_T * depth_m * depth_m * depth_m / 2.0) {
        if (!std::isfinite(colatitude_rad) || colatitude_rad <= 0.0 ||
            colatitude_rad >= std::numbers::pi / 2.0 ||
            !std::isfinite(depth_m) || depth_m <= 0.0 ||
            depth_m >= C::earthRadius_m || !std::isfinite(surfacePeak_T))
            throw std::invalid_argument("Invalid localized dipole parameters");
    }

    Vector3 anomaly(const Vector3& position) const {
        const Vector3 displacement = position - source;
        const double distance2 = displacement.norm2();
        if (!position.isFinite() || !std::isfinite(distance2) || distance2 <= 0.0)
            throw std::invalid_argument("Localized dipole undefined at source");
        const double distance = std::sqrt(distance2);
        return (displacement * (3.0 * axis.dot(displacement) / distance2) - axis) *
               (coefficient_T_m3 / (distance2 * distance));
    }

    Vector3 evaluate(const Vector3& position) const override {
        return background.evaluate(position) + anomaly(position);
    }
};

#endif // AURORA_LOCALIZEDDIPOLEFIELD_HPP
