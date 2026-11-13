#pragma once

namespace aurora::constants {
inline constexpr double speed_of_light_m_s = 299792458.0;
inline constexpr double elementary_charge_C = 1.602176634e-19;
inline constexpr double electron_charge_C = -elementary_charge_C;
inline constexpr double earth_radius_m = 6.371e6;
inline constexpr double electron_mass_kg = 9.1093837139e-31;
inline constexpr double electron_rest_energy_J =
    electron_mass_kg * speed_of_light_m_s * speed_of_light_m_s;

inline constexpr double eV_to_J(double energy_eV) {
    return energy_eV * elementary_charge_C;
}

inline constexpr double J_to_eV(double energy_J) {
    return energy_J / elementary_charge_C;
}
} // namespace aurora::constants
