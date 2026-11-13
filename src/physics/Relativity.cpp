#include "aurora/physics/Relativity.hpp"
#include "aurora/physics/Constants.hpp"

namespace aurora {
namespace {
double scaled_momentum(Vec3 p) {
    if (!finite(p))
        throw std::invalid_argument("Momentum must be finite");

    const double u = norm(p) / (constants::electron_mass_kg * constants::speed_of_light_m_s);
    if (!std::isfinite(u))
        throw std::overflow_error("Momentum is outside the supported double range");

    return u;
}
}

double lorentz_factor(Vec3 p) {
    return std::hypot(1.0, scaled_momentum(p));
}

Vec3 velocity_from_momentum(Vec3 p) {
    return (p / lorentz_factor(p)) / constants::electron_mass_kg;
}

double kinetic_energy_J(Vec3 p) {
    const double u = scaled_momentum(p);
    // gamma-1 = u * (u/(gamma+1))
    return constants::electron_rest_energy_J * u * (u / (std::hypot(1.0, u) + 1.0));
}

Vec3 momentum_from_velocity(Vec3 v) {
    if (!finite(v))
        throw std::invalid_argument("Velocity must be finite");

    const double beta = norm(v) / constants::speed_of_light_m_s;
    if (!(beta < 1.0))
        throw std::invalid_argument("Electron speed must be below c");

    const double gamma = 1.0 / std::sqrt((1.0 - beta) * (1.0 + beta));

    return v * (constants::electron_mass_kg * gamma);
}

Vec3 momentum_from_energy(double energy, Vec3 direction) {
    if (!std::isfinite(energy) || energy < 0.0)
        throw std::invalid_argument("Kinetic energy must be finite and non-negative");

    const Vec3 unit = normalized(direction);
    const double k = energy / constants::electron_rest_energy_J;
    const double magnitude = (constants::electron_mass_kg * constants::speed_of_light_m_s) *
        std::sqrt(k) * std::sqrt(k + 2.0);
    if (!std::isfinite(magnitude))
        throw std::overflow_error("Kinetic energy is outside the supported range");

    return unit * magnitude;
}
} // namespace aurora
