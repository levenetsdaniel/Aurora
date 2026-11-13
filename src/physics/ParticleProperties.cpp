#include "aurora/physics/ParticleProperties.hpp"
#include "aurora/physics/Constants.hpp"
#include <algorithm>

namespace aurora {
ParticleProperties particle_properties(const ParticleState& state, const MagneticField& field) {
    validate_state(state);
    const Vec3 b=field.at(state.position_m,state.time_s);
    const Vec3 axis=normalized(b); // Pitch and mu require a nonzero field.
    ParticleProperties d;
    d.field_T=norm(b);
    d.parallel_momentum_kg_m_s=dot(state.momentum_kg_m_s,axis);
    d.perpendicular_momentum_kg_m_s=norm(cross(state.momentum_kg_m_s,axis));
    const double p=norm(state.momentum_kg_m_s);
    if (p>0) {
        d.parallel_fraction=std::clamp(d.parallel_momentum_kg_m_s/p,-1.0,1.0);
        d.pitch_angle_rad=std::atan2(d.perpendicular_momentum_kg_m_s,d.parallel_momentum_kg_m_s);
    }
    d.magnetic_moment_J_T=(d.perpendicular_momentum_kg_m_s/constants::electron_mass_kg)
        *(d.perpendicular_momentum_kg_m_s/d.field_T)*0.5;
    if (!std::isfinite(d.magnetic_moment_J_T)) throw std::overflow_error("Magnetic moment overflow");
    return d;
}
double predicted_mirror_field_T(const ParticleState& state, const MagneticField& field) {
    const auto d=particle_properties(state,field);
    const double p=norm(state.momentum_kg_m_s);
    if (!(d.perpendicular_momentum_kg_m_s>0))
        throw std::domain_error("Mirror prediction requires nonzero perpendicular momentum");
    const double ratio=p/d.perpendicular_momentum_kg_m_s;
    const double result=d.field_T*ratio*ratio;
    if (!std::isfinite(result)) throw std::overflow_error("Mirror field prediction overflow");
    return result;
}
} // namespace aurora
