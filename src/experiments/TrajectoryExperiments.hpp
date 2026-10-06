#ifndef AURORA_TRAJECTORYEXPERIMENTS_HPP
#define AURORA_TRAJECTORYEXPERIMENTS_HPP

#include "../numerics/BorisIntegrator.hpp"
#include "../physics/Constants.hpp"
#include "../physics/Diagnostics.hpp"
#include "../physics/DipoleField.hpp"
#include "../physics/UniformField.hpp"
#include "../simulation/SimulationEngine.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <stdexcept>
#include <utility>
#include <vector>

namespace TrajectoryExperiments {

inline constexpr double field_T = 5.0e-5;
inline constexpr double energy_eV = 10000.0;

inline double momentumMagnitude() {
    return Relativity::momentumMagnitudeFromKineticEnergy_J(
        energy_eV * C::elementaryCharge_C, C::electronMass_kg);
}

inline double gyroPeriod(double magneticMagnitude_T) {
    const Vector3 momentum(momentumMagnitude(), 0.0, 0.0);
    return 2.0 * std::numbers::pi *
        Relativity::gammaFromMomentum(momentum, C::electronMass_kg) *
        C::electronMass_kg / (C::elementaryCharge_C * magneticMagnitude_T);
}

struct OrbitRow {
    int stepsPerPeriod;
    double dt_s;
    double exactPeriod_s;
    double measuredPeriod_s;
    double exactRadius_m;
    double radiusError_m;
    double positionError_m;
    double maximumRelativeEnergyError;
};

inline OrbitRow uniformOrbit(int stepsPerPeriod) {
    if (stepsPerPeriod < 4) throw std::invalid_argument("stepsPerPeriod must be >= 4");
    const UniformField field(Vector3(0.0, 0.0, field_T));
    const BorisIntegrator integrator(field, -C::elementaryCharge_C, C::electronMass_kg);
    constexpr double perpendicularFraction = 0.8;
    constexpr int periods = 4;
    const double momentum = momentumMagnitude();
    const double pPerpendicular = perpendicularFraction * momentum;
    const double pParallel = std::sqrt(1.0 - perpendicularFraction * perpendicularFraction) * momentum;
    const double gamma = Relativity::gammaFromMomentum(
        Vector3(momentum, 0.0, 0.0), C::electronMass_kg);
    const double omega = -C::elementaryCharge_C * field_T / (gamma * C::electronMass_kg);
    const double exactPeriod_s = gyroPeriod(field_T);
    const double dt_s = exactPeriod_s / stepsPerPeriod;
    const double vPerpendicular = pPerpendicular / (gamma * C::electronMass_kg);
    const double vParallel = pParallel / (gamma * C::electronMass_kg);
    const double exactRadius_m = std::abs(vPerpendicular / omega);
    const Vector3 center(0.0, -vPerpendicular / omega, 0.0);

    ParticleState state{Vector3{}, Vector3(pPerpendicular, 0.0, pParallel), 0.0};
    double totalPhase = 0.0;
    double maximumRelativeEnergyError = 0.0;
    double radiusError_m = 0.0;
    for (int i = 0; i < periods * stepsPerPeriod; ++i) {
        const ParticleState next = integrator.step(state, dt_s);
        const double crossZ = state.p.x * next.p.y - state.p.y * next.p.x;
        const double dotXY = state.p.x * next.p.x + state.p.y * next.p.y;
        totalPhase += std::atan2(crossZ, dotXY);
        state = next;
        const double radius = std::hypot(state.r.x - center.x, state.r.y - center.y);
        radiusError_m = std::max(radiusError_m, std::abs(radius - exactRadius_m));
        const double actualEnergy = Relativity::kineticEnergyFromMomentum_J(
            state.p, C::electronMass_kg) / C::elementaryCharge_C;
        maximumRelativeEnergyError = std::max(maximumRelativeEnergyError,
            std::abs(actualEnergy / energy_eV - 1.0));
    }
    const double time_s = periods * stepsPerPeriod * dt_s;
    const Vector3 exactPosition(
        vPerpendicular * std::sin(omega * time_s) / omega,
        vPerpendicular * (std::cos(omega * time_s) - 1.0) / omega,
        vParallel * time_s);
    return {stepsPerPeriod, dt_s, exactPeriod_s,
            2.0 * std::numbers::pi * time_s / std::abs(totalPhase),
            exactRadius_m, radiusError_m, (state.r - exactPosition).norm(),
            maximumRelativeEnergyError};
}

struct EventRow {
    int stepsPerPeriod;
    double dt_s;
    double exactTime_s;
    double simulatedTime_s;
    double timeError_s;
    double sphereResidual_m;
    std::uint64_t steps;
};

inline EventRow atmosphericCrossing(int stepsPerPeriod) {
    if (stepsPerPeriod < 4) throw std::invalid_argument("stepsPerPeriod must be >= 4");
    const double atmosphereRadius = C::earthRadius_m + 150000.0;
    const double initialHeightAboveBoundary_m = 2.0;
    const double momentum = momentumMagnitude();
    const double gamma = Relativity::gammaFromMomentum(
        Vector3(momentum, 0.0, 0.0), C::electronMass_kg);
    const double pPerpendicular = -0.8 * momentum;
    const double pParallel = 0.6 * momentum;
    const double vPerpendicular = pPerpendicular / (gamma * C::electronMass_kg);
    const double vParallel = pParallel / (gamma * C::electronMass_kg);
    const double omega = -C::elementaryCharge_C * field_T / (gamma * C::electronMass_kg);
    const auto exactRadius = [&](double time_s) {
        const double x = atmosphereRadius + initialHeightAboveBoundary_m +
                         vPerpendicular * std::sin(omega * time_s) / omega;
        const double y = vPerpendicular * (std::cos(omega * time_s) - 1.0) / omega;
        const double z = vParallel * time_s;
        return std::hypot(std::hypot(x, y), z);
    };
    double low_s = 0.0;
    double high_s = gyroPeriod(field_T) / 4.0;
    for (int i = 0; i < 80; ++i) {
        const double middle_s = 0.5 * (low_s + high_s);
        if (exactRadius(middle_s) > atmosphereRadius) low_s = middle_s;
        else high_s = middle_s;
    }
    const double exactTime_s = 0.5 * (low_s + high_s);

    const UniformField field(Vector3(0.0, 0.0, field_T));
    const double dt_s = gyroPeriod(field_T) / stepsPerPeriod;
    const SimulationEngine engine(field, -C::elementaryCharge_C, C::electronMass_kg,
        {dt_s, 150000.0, atmosphereRadius + 100.0, 2.0 * exactTime_s, 100000});
    const ParticleState initial{Vector3(atmosphereRadius + initialHeightAboveBoundary_m,
        0.0, 0.0),
        Vector3(pPerpendicular, 0.0, pParallel), 0.0};
    const SimulationResult result = engine.run(initial);
    if (result.outcome != ParticleOutcome::Precipitated)
        throw std::runtime_error("control particle did not cross the atmospheric sphere");
    return {stepsPerPeriod, dt_s, exactTime_s, result.eventState.time_s,
            std::abs(result.eventState.time_s - exactTime_s),
            std::abs(result.eventState.r.norm() - atmosphereRadius), result.steps};
}

struct DiagnosticRow {
    double time_s;
    double energy_eV;
    double relativeEnergyError;
    double magneticMagnitude_T;
    double parallelMomentum_kg_m_s;
    double magneticMoment_J_T;
};

inline DiagnosticRow diagnosticRow(const ParticleState& state,
                                   const MagneticField& field) {
    const Vector3 magnetic = field.evaluate(state.r);
    const auto diagnostics = evaluateDiagnostics(state.p, magnetic, C::electronMass_kg);
    return {state.time_s, diagnostics.energy_eV,
            diagnostics.energy_eV / energy_eV - 1.0, magnetic.norm(),
            diagnostics.parallelMomentum_kg_m_s, diagnostics.magneticMoment_J_T};
}

inline std::vector<DiagnosticRow> energyTrace(const MagneticField& field,
                                               const ParticleState& initial,
                                               double dt_s, int steps,
                                               int stride) {
    if (steps <= 0 || stride <= 0) throw std::invalid_argument("invalid trace length");
    const BorisIntegrator integrator(field, -C::elementaryCharge_C, C::electronMass_kg);
    std::vector<DiagnosticRow> rows;
    ParticleState state = initial;
    rows.push_back(diagnosticRow(state, field));
    for (int i = 1; i <= steps; ++i) {
        state = integrator.step(state, dt_s);
        if (i % stride == 0 || i == steps) rows.push_back(diagnosticRow(state, field));
    }
    return rows;
}

struct MirrorExperiment {
    ParticleOutcome outcome;
    double initialMagneticMagnitude_T;
    double mirrorMagneticMagnitude_T;
    double predictedMirrorMagneticMagnitude_T;
    double relativeEnergyError;
    std::uint64_t steps;
    std::vector<DiagnosticRow> trace;
};

inline MirrorExperiment mirrorTrajectory(int stepsPerPeriod = 50,
                                         int traceStride = 100) {
    if (stepsPerPeriod < 4 || traceStride < 1)
        throw std::invalid_argument("invalid mirror resolution");
    const DipoleField field;
    constexpr double colatitude = 0.3;
    constexpr double pitch = 80.0 * std::numbers::pi / 180.0;
    const double radius = C::earthRadius_m + 450000.0;
    const Vector3 position(radius * std::sin(colatitude), 0.0,
                           radius * std::cos(colatitude));
    const Vector3 magnetic = field.evaluate(position);
    const Vector3 along = magnetic * (-1.0 / magnetic.norm());
    const Vector3 across = Vector3(0.0, 1.0, 0.0);
    const Vector3 momentum = (along * std::cos(pitch) + across * std::sin(pitch)) *
                             momentumMagnitude();
    const ParticleState initial{position, momentum, 0.0};
    const double dt_s = gyroPeriod(magnetic.norm()) / stepsPerPeriod;
    const SimulationEngine engine(field, -C::elementaryCharge_C, C::electronMass_kg,
        {dt_s, 150000.0, C::earthRadius_m + 1000000.0, 0.03, 2000000});
    std::vector<DiagnosticRow> trace;
    std::uint64_t samples = 0;
    const SimulationResult result = engine.run(initial, [&](const ParticleState& state) {
        if (samples++ % static_cast<std::uint64_t>(traceStride) == 0)
            trace.push_back(diagnosticRow(state, field));
    });
    trace.push_back(diagnosticRow(result.eventState, field));
    const double initialB = magnetic.norm();
    const double eventB = field.evaluate(result.eventState.r).norm();
    const double finalEnergy = Relativity::kineticEnergyFromMomentum_J(
        result.eventState.p, C::electronMass_kg) / C::elementaryCharge_C;
    return {result.outcome, initialB, eventB, initialB / std::pow(std::sin(pitch), 2),
            finalEnergy / energy_eV - 1.0, result.steps, std::move(trace)};
}

} // namespace TrajectoryExperiments

#endif // AURORA_TRAJECTORYEXPERIMENTS_HPP
