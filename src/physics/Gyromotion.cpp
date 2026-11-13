#include "aurora/physics/Gyromotion.hpp"
#include "aurora/physics/Constants.hpp"
#include "aurora/physics/Relativity.hpp"
#include <numbers>

namespace aurora {
double gyrofrequency_rad_s(Vec3 p, Vec3 b) {
    if (!finite(b) || !std::isfinite(norm(b))) throw std::invalid_argument("Invalid field");
    const double omega = (constants::elementary_charge_C / constants::electron_mass_kg)
        * (norm(b) / lorentz_factor(p));
    if (!std::isfinite(omega)) throw std::overflow_error("Gyrofrequency overflow");
    return omega;
}
double gyroperiod_s(Vec3 p, Vec3 b) {
    const double omega = gyrofrequency_rad_s(p,b);
    if (!(omega > 0.0)) throw std::invalid_argument("Gyroperiod is undefined for B=0");
    return 2.0*std::numbers::pi/omega;
}
double larmor_radius_m(Vec3 p, Vec3 b) {
    const Vec3 axis = normalized(b);
    return norm(p-axis*dot(p,axis)) / (constants::elementary_charge_C * norm(b));
}
} // namespace aurora
