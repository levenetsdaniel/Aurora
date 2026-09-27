#include "../src/ensemble/ElectronSource.hpp"
#include "../src/ensemble/PrecipitationCollector.hpp"
#include "../src/physics/UniformField.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

template <typename Action>
void requireInvalid(Action action, const char* message) {
    try {
        action();
    } catch (const std::invalid_argument&) {
        return;
    }
    throw std::runtime_error(message);
}

void checkEventStateConversion() {
    const UniformField field(Vector3(-1e-9, 0.0, 0.0));
    const double radius_m = C::earthRadius_m + 150000.0;
    const double momentum = Relativity::momentumMagnitudeFromKineticEnergy_J(
        2000.0 * C::elementaryCharge_C, C::electronMass_kg);
    const ParticleState event{Vector3(radius_m, 0.0, 0.0),
                              Vector3(-std::cos(0.3), std::sin(0.3), 0.0) * momentum,
                              0.002};
    const InitialCondition initial{37,
        {Vector3(radius_m + 50000.0, 0.0, 0.0), Vector3(-momentum, 0.0, 0.0), 0.0},
        1000.0, 0.0, 0.0, 2.5};
    const EnsembleRecord record{37, initial,
        {ParticleOutcome::Precipitated, event, 12, 0.0, 0.0, 0.0}};

    const PrecipitationCollector collector;
    const auto collected = collector.collect(record, field, C::electronMass_kg);
    require(collected.has_value(), "precipitated particle was discarded");
    require(collected->id == 37 && collected->weight == 2.5 &&
            collected->state.r.x == event.r.x && collected->state.p.x == event.p.x &&
            collected->state.time_s == event.time_s,
            "collector did not preserve id, event state, or marker weight");
    require(std::abs(collected->energy_eV - 2000.0) < 1e-9 &&
            collected->energy_eV != initial.energy_eV,
            "collector used initial metadata instead of event momentum for energy");
    require(std::abs(collected->incidenceCosine - std::cos(0.3)) < 1e-12,
            "collector computed the wrong inward incidence cosine");

    for (ParticleOutcome outcome : {ParticleOutcome::Mirrored,
                                    ParticleOutcome::Escaped,
                                    ParticleOutcome::Unresolved}) {
        EnsembleRecord other = record;
        other.result.outcome = outcome;
        require(!collector.collect(other, field, C::electronMass_kg).has_value(),
                "non-precipitated particle entered atmospheric input");
    }

    EnsembleRecord invalid = record;
    invalid.initial.weight = -1.0;
    requireInvalid([&] { (void)collector.collect(invalid, field, C::electronMass_kg); },
                   "negative marker weight was accepted");
    invalid = record;
    invalid.result.eventState.p = Vector3(momentum, 0.0, 0.0);
    requireInvalid([&] { (void)collector.collect(invalid, field, C::electronMass_kg); },
                   "outgoing event momentum was accepted");
    invalid = record;
    invalid.result.eventState.p = Vector3();
    requireInvalid([&] { (void)collector.collect(invalid, field, C::electronMass_kg); },
                   "zero event speed was accepted");
    requireInvalid([&] { (void)collector.collect(record, field, 0.0); },
                   "zero particle mass was accepted");
}

void checkSourceToAtmosphereBoundary() {
    const UniformField field(Vector3(-1e-9, 0.0, 0.0));
    const Vector3 start(C::earthRadius_m + 200000.0, 0.0, 0.0);
    const IsotropicElectronSource source(field, start, 0.0, C::electronMass_kg,
                                         1000.0, 2000.0, 1.0, 1.0, 100, true);
    const SimulationConfig config{1e-5, 150000.0,
        C::earthRadius_m + 300000.0, 0.01, 1000};
    const SimulationEngine engine(field, 0.0, C::electronMass_kg, config);
    const EnsembleRunner runner(engine);
    const PrecipitationCollector collector;
    Random rng(123);
    InitialCondition incoming;
    bool found = false;
    for (std::uint64_t id = 0; id < 1000; ++id) {
        const InitialCondition candidate = source.sample(rng, id);
        if (candidate.pitch_rad < 0.2) {
            incoming = candidate;
            found = true;
            break;
        }
    }
    require(found, "test source did not sample a near-radial electron");

    const auto records = runner.run({incoming});
    require(records.size() == 1 &&
            records.front().result.outcome == ParticleOutcome::Precipitated,
            "inward source electron did not reach the atmosphere");
    const auto collected = collector.collect(records.front(), field, C::electronMass_kg);
    require(collected.has_value(), "precipitated source electron was not collected");
    require(std::abs(collected->state.r.norm() -
                     (C::earthRadius_m + config.atmosphereHeight_m)) < 1e-3,
            "atmospheric entry state was not localized to the boundary");
    require(collected->weight == incoming.weight && collected->id == incoming.id &&
            collected->incidenceCosine > 0.0,
            "source metadata or inward incidence was lost");
}

} // namespace

int main() {
    try {
        checkEventStateConversion();
        checkSourceToAtmosphereBoundary();
        std::cout << "PrecipitationCollector checks passed\n";
    } catch (const std::exception& error) {
        std::cerr << "PrecipitationCollector check failed: " << error.what() << '\n';
        return 1;
    }
}
