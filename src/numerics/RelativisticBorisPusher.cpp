#include "aurora/numerics/RelativisticBorisPusher.hpp"
#include "aurora/physics/Constants.hpp"
#include "aurora/physics/Relativity.hpp"
#include <limits>

namespace aurora {
RelativisticBorisPusher::RelativisticBorisPusher(double dt_s) : dt_s_(dt_s) {
    if (!std::isfinite(dt_s) || !(dt_s > 0.0))
        throw std::invalid_argument("Time step must be finite and positive");
}

void RelativisticBorisPusher::check_step(const State& s) const {
    if (s.dt_s_ != dt_s_)
        throw std::invalid_argument("Cannot change the step of a staggered trajectory");
}

Vec3 RelativisticBorisPusher::rotate(Vec3 p, Vec3 b, double interval) {
    if (!finite(b))
        throw std::runtime_error("Field returned a non-finite value");

    const double gamma = lorentz_factor(p);
    // t_vec = q_e B * interval / (2 m_e gamma), безразмерный вектор вращения.
    const Vec3 t = b *
        ((constants::electron_charge_C / constants::electron_mass_kg) * (0.5 * interval / gamma));
    const double t2 = dot(t, t);
    if (!finite(t) || !std::isfinite(t2))
        throw std::overflow_error("Boris rotation overflow");

    // s_vec = 2 t_vec / (1 + |t_vec|^2).
    const Vec3 s = t * (2.0 / (1.0 + t2));

    // p_prime = p + p × t_vec; p_new = p + p_prime × s_vec.
    const Vec3 result = p + cross(p + cross(p, t), s);
    if (!finite(result))
        throw std::overflow_error("Non-finite momentum after Boris rotation");

    return result;
}

RelativisticBorisPusher::State RelativisticBorisPusher::initialize(
    const ParticleState& initial, const MagneticField& field) const {
    validate_state(initial);

    State s;
    s.position_m_ = initial.position_m;
    // p^(-1/2) = R(B(r^0,t_0), -dt/2) p^0: подготовка импульса на половинном слое.
    s.momentum_before_ =
        rotate(initial.momentum_kg_m_s, field.at(initial.position_m, initial.time_s), -0.5 * dt_s_);
    s.initial_time_s_ = initial.time_s;
    s.dt_s_ = dt_s_;

    return s;
}

void RelativisticBorisPusher::advance(State& s, const MagneticField& field) const {
    check_step(s);
    if (s.steps_ == std::numeric_limits<std::size_t>::max())
        throw std::overflow_error("Step counter overflow");

    // p^(n+1/2) = R(B(r^n,t_n), dt) p^(n-1/2).
    const Vec3 p_after = rotate(s.momentum_before_, field.at(s.position_m_, s.time_s()), dt_s_);

    // v^(n+1/2) = p^(n+1/2)/(gamma m_e); r^(n+1) = r^n + dt * v^(n+1/2).
    const Vec3 position = s.position_m_ + velocity_from_momentum(p_after) * dt_s_;
    const double next_time = s.initial_time_s_ + static_cast<double>(s.steps_ + 1) * dt_s_;
    if (!finite(position) || !std::isfinite(next_time) || !(next_time > s.time_s()))
        throw std::overflow_error("Time or position cannot advance in double precision");

    s.position_m_ = position;
    s.momentum_before_ = p_after;
    ++s.steps_;
}

ParticleState RelativisticBorisPusher::synchronize(
    const State& s, const MagneticField& field) const {
    check_step(s);

    // p_sync^n = R(B(r^n,t_n), dt/2) p^(n-1/2).
    return {s.position_m_,
        rotate(s.momentum_before_, field.at(s.position_m_, s.time_s()), 0.5 * dt_s_), s.time_s()};
}
} // namespace aurora
