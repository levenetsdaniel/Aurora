#ifndef AURORA_ELECTRONSCENEDATA_HPP
#define AURORA_ELECTRONSCENEDATA_HPP

#include "../physics/Constants.hpp"
#include "../physics/MagneticField.hpp"
#include "../ensemble/InitialCondition.hpp"
#include "../simulation/SimulationEngine.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

struct ElectronPath {
    std::uint64_t id{};
    double energy_eV{};
    double weight_s1{};
    ParticleOutcome outcome{ParticleOutcome::Unresolved};
    ParticleState terminalState;
    std::vector<Vector3> positions_m;
    double initialPitch_rad{};
    double initialGyroPhase_rad{};
    ParticleState initialState;
};

struct ElectronSceneData {
    double atmosphereHeight_m{};
    std::vector<ElectronPath> paths;
    std::vector<std::vector<Vector3>> magneticLines_m;

    std::size_t precipitationCount() const {
        return static_cast<std::size_t>(std::count_if(paths.begin(), paths.end(),
            [](const ElectronPath& path) {
                return path.outcome == ParticleOutcome::Precipitated;
            }));
    }

    std::size_t countOutcome(ParticleOutcome outcome) const {
        return static_cast<std::size_t>(std::count_if(paths.begin(), paths.end(),
            [outcome](const ElectronPath& path) { return path.outcome == outcome; }));
    }
};

class ElectronSceneSampler {
private:
    const SimulationEngine& engine;
    const MagneticField& field;

    static Vector3 fieldDirection(const MagneticField& field,
                                  const Vector3& position_m, double sign) {
        const Vector3 value = field.evaluate(position_m);
        const double magnitude = value.norm();
        if (!value.isFinite() || !std::isfinite(magnitude) || magnitude <= 0.0) {
            throw std::runtime_error("ElectronSceneSampler: invalid magnetic field");
        }
        return value * (sign / magnitude);
    }

    static std::vector<Vector3> traceHalf(const MagneticField& field,
                                           Vector3 seed_m, double sign,
                                           double step_m, std::size_t maxSteps) {
        std::vector<Vector3> points{seed_m};
        for (std::size_t index = 0; index < maxSteps; ++index) {
            const Vector3 k1 = fieldDirection(field, seed_m, sign);
            const Vector3 k2 = fieldDirection(field, seed_m + k1 * (step_m / 2), sign);
            const Vector3 k3 = fieldDirection(field, seed_m + k2 * (step_m / 2), sign);
            const Vector3 k4 = fieldDirection(field, seed_m + k3 * step_m, sign);
            const Vector3 next = seed_m +
                (k1 + k2 * 2 + k3 * 2 + k4) * (step_m / 6);
            const double radius = next.norm();
            if (!next.isFinite() || !std::isfinite(radius) ||
                radius <= C::earthRadius_m || radius >= 2.0 * C::earthRadius_m) {
                break;
            }
            points.push_back(next);
            seed_m = next;
        }
        return points;
    }

public:
    ElectronSceneSampler(const SimulationEngine& engine, const MagneticField& field)
        : engine(engine), field(field) {}

    ElectronPath traceElectron(const InitialCondition& initial,
                               std::size_t sampleStride = 250) const {
        if (sampleStride == 0 || !std::isfinite(initial.energy_eV) ||
            initial.energy_eV <= 0.0 || !std::isfinite(initial.weight) ||
            initial.weight < 0.0) {
            throw std::invalid_argument("ElectronSceneSampler: invalid marker or stride");
        }
        ElectronPath path{initial.id, initial.energy_eV, initial.weight};
        path.initialPitch_rad = initial.pitch_rad;
        path.initialGyroPhase_rad = initial.gyroPhase_rad;
        path.initialState = initial.state;
        std::uint64_t seen = 0;
        const SimulationResult result = engine.run(initial.state,
            [&](const ParticleState& state) {
                if (seen % sampleStride == 0) path.positions_m.push_back(state.r);
                ++seen;
            });
        path.outcome = result.outcome;
        path.terminalState = result.eventState;
        const Vector3 terminal = result.eventState.r;
        if (path.positions_m.empty() ||
            (path.positions_m.back() - terminal).norm() > 0.0) {
            path.positions_m.push_back(terminal);
        }
        return path;
    }

    std::vector<Vector3> traceMagneticLine(const Vector3& seed_m,
                                            double step_m = 20000.0,
                                            std::size_t maxSteps = 1500) const {
        if (!seed_m.isFinite() || seed_m.norm() <= C::earthRadius_m ||
            !std::isfinite(step_m) || step_m <= 0.0 || maxSteps == 0) {
            throw std::invalid_argument("ElectronSceneSampler: invalid field line seed");
        }
        auto backward = traceHalf(field, seed_m, -1.0, step_m, maxSteps);
        auto forward = traceHalf(field, seed_m, 1.0, step_m, maxSteps);
        std::reverse(backward.begin(), backward.end());
        backward.insert(backward.end(), forward.begin() + 1, forward.end());
        return backward;
    }
};

#endif // AURORA_ELECTRONSCENEDATA_HPP
