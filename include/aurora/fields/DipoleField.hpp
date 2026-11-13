#pragma once

#include "aurora/fields/MagneticField.hpp"
#include "aurora/physics/Constants.hpp"

namespace aurora {
// Axis is the direction of the magnetic dipole MOMENT (not the north magnetic pole).
class DipoleField final : public MagneticField {
public:
    explicit DipoleField(double equatorial_surface_T = 3.12e-5,
                         Vec3 moment_axis = {0,0,-1},
                         double reference_radius_m = constants::earth_radius_m,
                         double minimum_radius_m = 0.5*constants::earth_radius_m);
    Vec3 at(Vec3 position_m, double time_s) const override;
    Vec3 magnitude_gradient(Vec3 position_m) const; // grad |B|, T/m
    Vec3 axis() const { return axis_; }
    double reference_radius_m() const { return radius_m_; }
    double minimum_radius_m() const { return minimum_m_; }
private:
    double equatorial_T_, radius_m_, minimum_m_;
    Vec3 axis_;
    double checked_radius(Vec3 position_m) const;
};
} // namespace aurora
