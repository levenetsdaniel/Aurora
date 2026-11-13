#include "reference/UniformFieldSolution.hpp"
#include "aurora/physics/Gyromotion.hpp"
#include "aurora/physics/Constants.hpp"
#include "aurora/physics/Relativity.hpp"

namespace aurora::test_reference {
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
} // namespace aurora::test_reference
