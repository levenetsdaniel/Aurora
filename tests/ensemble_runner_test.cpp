#include "../src/ensemble/EnsembleRunner.hpp"
#include "../src/physics/Constants.hpp"
#include "../src/physics/Relativity.hpp"
#include "../src/physics/UniformField.hpp"

#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

bool sameVector(const Vector3& left, const Vector3& right) {
    return left.x == right.x && left.y == right.y && left.z == right.z;
}

bool sameState(const ParticleState& left, const ParticleState& right) {
    return sameVector(left.r, right.r) && sameVector(left.p, right.p) &&
           left.time_s == right.time_s;
}

bool sameResult(const SimulationResult& left, const SimulationResult& right) {
    return left.outcome == right.outcome && sameState(left.eventState, right.eventState) &&
           left.steps == right.steps && left.initialEnergy_eV == right.initialEnergy_eV &&
           left.initialPitch_deg == right.initialPitch_deg &&
           left.maxRelativeEnergyError == right.maxRelativeEnergyError;
}

std::vector<InitialCondition> makeGrid() {
    const Vector3 start(C::earthRadius_m + 200000.0, 0.0, 0.0);
    std::vector<InitialCondition> conditions;
    std::uint64_t id = 100;

    for (double energy_eV : {1000.0, 2000.0}) {
        const double momentum = Relativity::momentumMagnitudeFromKineticEnergy_J(
            energy_eV * C::elementaryCharge_C, C::electronMass_kg);
        for (double pitch_rad : {0.2, 0.4}) {
            for (double phase_rad : {0.0, std::numbers::pi}) {
                const Vector3 direction(-std::cos(pitch_rad),
                                        std::sin(pitch_rad) * std::cos(phase_rad),
                                        std::sin(pitch_rad) * std::sin(phase_rad));
                conditions.push_back({id, {start, direction * momentum, 0.0},
                                      energy_eV, pitch_rad, phase_rad, 0.5});
                id += 7;
            }
        }
    }

    const double momentum = Relativity::momentumMagnitudeFromKineticEnergy_J(
        1000.0 * C::elementaryCharge_C, C::electronMass_kg);
    conditions.push_back({id, {start, Vector3(momentum, 0.0, 0.0), 0.0},
                          1000.0, 0.0, 0.0, 2.0});
    conditions.push_back({id + 7, {start, Vector3(), 0.0},
                          0.0, 0.0, 0.0, 0.0});
    return conditions;
}

void checkEnsemble() {
    const UniformField field(Vector3(1e-9, 0.0, 0.0));
    const SimulationConfig config{1e-5, 150000.0,
        C::earthRadius_m + 300000.0, 0.01, 1000};
    const SimulationEngine engine(field, 0.0, C::electronMass_kg, config);
    const EnsembleRunner runner(engine);
    const std::vector<InitialCondition> conditions = makeGrid();

    require(runner.run({}).empty(), "empty ensemble did not produce an empty result");

    const std::vector<EnsembleRecord> records = runner.run(conditions);
    const std::vector<EnsembleRecord> repeated = runner.run(conditions);
    require(records.size() == conditions.size() && repeated.size() == records.size(),
            "ensemble lost or added particles");

    for (std::size_t index = 0; index < records.size(); ++index) {
        const InitialCondition& initial = conditions[index];
        const EnsembleRecord& record = records[index];
        require(record.id == initial.id && record.initial.id == initial.id &&
                sameState(record.initial.state, initial.state) &&
                record.initial.energy_eV == initial.energy_eV &&
                record.initial.pitch_rad == initial.pitch_rad &&
                record.initial.gyroPhase_rad == initial.gyroPhase_rad &&
                record.initial.weight == initial.weight,
                "ensemble record lost initial-condition metadata or order");
        require(sameResult(record.result, engine.run(initial.state)),
                "ensemble result differs from an individual run");
        require(repeated[index].id == record.id &&
                sameResult(repeated[index].result, record.result),
                "repeated ensemble run was not reproducible");
    }

    for (std::size_t index = 0; index < 8; ++index) {
        require(records[index].result.outcome == ParticleOutcome::Precipitated,
                "inward grid particle did not precipitate");
    }
    require(records[8].result.outcome == ParticleOutcome::Escaped,
            "outward particle did not escape");
    require(records[9].result.outcome == ParticleOutcome::Unresolved,
            "stationary particle was not unresolved");

    std::vector<InitialCondition> invalid = conditions;
    invalid[3].state.r.x = C::earthRadius_m;
    try {
        (void)runner.run(invalid);
    } catch (const std::invalid_argument&) {
        return;
    }
    throw std::runtime_error("invalid initial state did not propagate an engine error");
}

} // namespace

int main() {
    try {
        checkEnsemble();
        std::cout << "EnsembleRunner checks passed\n";
    } catch (const std::exception& error) {
        std::cerr << "EnsembleRunner check failed: " << error.what() << '\n';
        return 1;
    }
}
