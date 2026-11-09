#include "aurora/diagnostics/UniformFieldSolution.hpp"
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
ParticleState uniform_field_solution(const ParticleState& initial, Vec3 b, double time) {
    validate_state(initial);
    if (!std::isfinite(time)) throw std::invalid_argument("Reference time must be finite");
    const double dt = time-initial.time_s;
    const Vec3 v = velocity_from_momentum(initial.momentum_kg_m_s);
    const double omega = gyrofrequency_rad_s(initial.momentum_kg_m_s,b);
    if (omega == 0.0) return {initial.position_m+v*dt, initial.momentum_kg_m_s, time};
    const Vec3 axis = normalized(b);
    const Vec3 parallel = axis*dot(initial.momentum_kg_m_s,axis);
    const Vec3 perpendicular = initial.momentum_kg_m_s-parallel;
    // q_e < 0: dp/dt = q p×B/(gamma*m) = +omega (axis×p).
    const double angle = omega*dt;
    const double sin_angle = std::sin(angle);
    // 1-cos(angle) without cancellation for small angles.
    const double one_minus_cos = 2.0*std::pow(std::sin(0.5*angle),2);
    const double gamma_mass = lorentz_factor(initial.momentum_kg_m_s)*constants::electron_mass_kg;
    const Vec3 p = parallel + perpendicular*std::cos(angle) + cross(axis,perpendicular)*sin_angle;
    const Vec3 displacement = (parallel*dt + perpendicular*(sin_angle/omega)
        + cross(axis,perpendicular)*(one_minus_cos/omega)) / gamma_mass;
    return {initial.position_m+displacement,p,time};
}
} // namespace aurora
