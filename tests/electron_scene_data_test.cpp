#include "../src/visualization/ElectronSceneData.hpp"
#include "../src/ensemble/ArcElectronSource.hpp"
#include "../src/physics/DipoleField.hpp"
#include "../src/physics/Relativity.hpp"
#include "../src/physics/UniformField.hpp"

#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void checkSceneData() {
    const UniformField field({1e-9, 0, 0});
    const SimulationConfig config{1e-5, 150000,
        C::earthRadius_m + 300000, 0.01, 1000};
    const SimulationEngine engine(field, -C::elementaryCharge_C,
                                  C::electronMass_kg, config);
    const ElectronSceneSampler sampler(engine, field);
    const double momentum = Relativity::momentumMagnitudeFromKineticEnergy_J(
        1000 * C::elementaryCharge_C, C::electronMass_kg);
    const InitialCondition initial{42,
        {{C::earthRadius_m + 200000, 0, 0}, {-momentum, 0, 0}, 0},
        1000, 0, 0, 10};
    const ElectronPath path = sampler.traceElectron(initial, 5);
    require(path.id == 42 && path.energy_eV == 1000 && path.weight_s1 == 10 &&
            path.outcome == ParticleOutcome::Precipitated &&
            path.positions_m.size() > 2 &&
            path.positions_m.front().x == initial.state.r.x &&
            std::abs(path.positions_m.back().norm() -
                     (C::earthRadius_m + 150000)) < 1e-3,
            "electron path lost metadata or event boundary");
    const auto line = sampler.traceMagneticLine(initial.state.r, 20000, 100);
    require(line.size() > 2 &&
            (line.front() - initial.state.r).norm() > 0 &&
            (line.back() - initial.state.r).norm() > 0,
            "magnetic field line did not trace in both directions");
    const ElectronSceneData scene{150000, {path}, {line}};
    require(scene.precipitationCount() == 1 &&
            scene.countOutcome(ParticleOutcome::Mirrored) == 0 &&
            path.initialPitch_rad == initial.pitch_rad &&
            path.initialGyroPhase_rad == initial.gyroPhase_rad &&
            path.initialState.r.x == initial.state.r.x &&
            path.initialState.p.x == initial.state.p.x,
            "scene precipitation count is wrong");
}

void checkWideAngleOutcomes() {
    const DipoleField field;
    const SimulationConfig config{2e-8, 150000,
        C::earthRadius_m + 800000, 0.05, 2500000};
    const SimulationEngine engine(field, -C::elementaryCharge_C,
                                  C::electronMass_kg, config);
    const ArcElectronSource source(field, true, 0.0, std::numbers::pi / 2.0,
                                    1e20, 32);
    Random random(20260928);
    ElectronSceneData scene;
    for (std::uint64_t id = 1; id <= 32; ++id) {
        scene.paths.push_back(ElectronSceneSampler(engine, field).traceElectron(
            source.sample(random, id), 250));
    }
    require(scene.countOutcome(ParticleOutcome::Precipitated) > 0 &&
            scene.countOutcome(ParticleOutcome::Mirrored) > 0 &&
            scene.countOutcome(ParticleOutcome::Unresolved) == 0 &&
            scene.countOutcome(ParticleOutcome::Precipitated) +
                scene.countOutcome(ParticleOutcome::Mirrored) +
                scene.countOutcome(ParticleOutcome::Escaped) +
                scene.countOutcome(ParticleOutcome::Unresolved) == scene.paths.size(),
            "wide-angle arc did not contain both precipitating and mirrored particles");
}

} // namespace

int main() {
    try {
        checkSceneData();
        checkWideAngleOutcomes();
        std::cout << "ElectronSceneData checks passed\n";
    } catch (const std::exception& error) {
        std::cerr << "ElectronSceneData check failed: " << error.what() << '\n';
        return 1;
    }
}
