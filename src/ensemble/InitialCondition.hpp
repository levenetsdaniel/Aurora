#ifndef AURORA_INITIALCONDITION_HPP
#define AURORA_INITIALCONDITION_HPP

#include "../simulation/Particle.hpp"

#include <cstdint>

struct InitialCondition {
    std::uint64_t id{};
    ParticleState state;
    double energy_eV{};
    double pitch_rad{};
    double gyroPhase_rad{};
    // Physical electrons per second represented by this marker in a fixed-size batch.
    double weight{1.0};
};

#endif // AURORA_INITIALCONDITION_HPP
