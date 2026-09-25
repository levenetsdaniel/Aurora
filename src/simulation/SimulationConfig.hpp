#ifndef AURORA_SIMULATIONCOMFIG_HPP
#define AURORA_SIMULATIONCOMFIG_HPP

#include <cstdint>

struct SimulationConfig {
    double dt_s;
    double atmosphereHeight_m;
    double escapeRadius_m;
    double maxTime_s;
    std::uint64_t maxSteps;
};


#endif //AURORA_SIMULATIONCOMFIG_HPP
