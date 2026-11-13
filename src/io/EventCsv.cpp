#include "aurora/io/EventCsv.hpp"
#include "aurora/physics/Constants.hpp"
#include "aurora/physics/Relativity.hpp"
#include <iomanip>
#include <ostream>

namespace aurora {
void write_event_header(std::ostream& out) {
    out << "type,time_s,x_m,y_m,z_m,B_T,bracket_width_s,altitude_m,kinetic_energy_J,"
           "px_kg_m_s,py_kg_m_s,pz_kg_m_s\n";

    if (!out)
        throw std::runtime_error("Cannot write event header");
}

void write_event_row(std::ostream& out, const TrajectoryEvent& event, const MagneticField& field,
    std::string_view inner_boundary_name) {
    const auto& s = event.state;
    const std::string_view name = event.type == EventType::Mirror
        ? "mirror"
        : (event.type == EventType::InnerBoundary ? inner_boundary_name : "outer_boundary");

    out << std::setprecision(17) << name << ',' << s.time_s << ',' << s.position_m.x << ','
        << s.position_m.y << ',' << s.position_m.z << ',' << norm(field.at(s.position_m, s.time_s))
        << ',' << event.bracket_width_s << ',' << norm(s.position_m) - constants::earth_radius_m
        << ',' << kinetic_energy_J(s.momentum_kg_m_s) << ',' << s.momentum_kg_m_s.x << ','
        << s.momentum_kg_m_s.y << ',' << s.momentum_kg_m_s.z << '\n';

    if (!out)
        throw std::runtime_error("Cannot write event row");
}
} // namespace aurora
