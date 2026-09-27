#ifndef AURORA_ENSEMBLERUNNER_HPP
#define AURORA_ENSEMBLERUNNER_HPP

#include "InitialCondition.hpp"
#include "../simulation/SimulationEngine.hpp"

#include <cstdint>
#include <vector>

struct EnsembleRecord {
    std::uint64_t id{};
    InitialCondition initial;
    SimulationResult result;
};

class EnsembleRunner {
private:
    const SimulationEngine& engine;

public:
    explicit EnsembleRunner(const SimulationEngine& engine) : engine(engine) {}

    std::vector<EnsembleRecord> run(
        const std::vector<InitialCondition>& initialConditions) const {
        std::vector<EnsembleRecord> records;
        records.reserve(initialConditions.size());

        for (const InitialCondition& initial : initialConditions) {
            records.push_back({initial.id, initial, engine.run(initial.state)});
        }

        return records;
    }
};

#endif // AURORA_ENSEMBLERUNNER_HPP
