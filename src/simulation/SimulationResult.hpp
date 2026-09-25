#ifndef AURORA_SIMULATIONDETECTOR_HPP
#define AURORA_SIMULATIONDETECTOR_HPP

#include <cstdint>

#include "Particle.hpp"


struct SimulationResult {
    ParticleOutcome outcome;
    ParticleState eventState;
    std::uint64_t steps;
    double initialEnergy_eV;
    double initialPitch_deg;
    double maxRelativeEnergyError;
};


#endif //AURORA_SIMULATIONDETECTOR_HPP
