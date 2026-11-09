#include "aurora/io/TrajectoryCsv.hpp"
#include "aurora/physics/Relativity.hpp"
#include <iomanip>
#include <ostream>

namespace aurora {
void write_trajectory_header(std::ostream& out) {
    out << "step,time_s,x_m,y_m,z_m,px_kg_m_s,py_kg_m_s,pz_kg_m_s,"
           "vx_m_s,vy_m_s,vz_m_s,gamma,kinetic_energy_J\n";
    if (!out) throw std::runtime_error("Cannot write trajectory header");
}
void write_trajectory_row(std::ostream& out, std::size_t step, const ParticleState& s) {
    const auto v = velocity_from_momentum(s.momentum_kg_m_s);
    out << std::setprecision(17) << step << ',' << s.time_s << ','
        << s.position_m.x << ',' << s.position_m.y << ',' << s.position_m.z << ','
        << s.momentum_kg_m_s.x << ',' << s.momentum_kg_m_s.y << ',' << s.momentum_kg_m_s.z << ','
        << v.x << ',' << v.y << ',' << v.z << ',' << lorentz_factor(s.momentum_kg_m_s) << ','
        << kinetic_energy_J(s.momentum_kg_m_s) << '\n';
    if (!out) throw std::runtime_error("Cannot write trajectory sample");
}
} // namespace aurora
