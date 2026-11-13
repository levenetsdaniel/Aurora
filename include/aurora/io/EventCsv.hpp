#pragma once

#include "aurora/simulation/TrajectoryRunner.hpp"
#include <iosfwd>
#include <string_view>

namespace aurora {
void write_event_header(std::ostream& out);
// The caller names the inner boundary according to the physical experiment.
void write_event_row(std::ostream& out, const TrajectoryEvent& event, const MagneticField& field,
    std::string_view inner_boundary_name = "inner_boundary");
} // namespace aurora
