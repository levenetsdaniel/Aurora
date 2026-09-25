#include "../src/numerics/BorisIntegrator.hpp"
#include "../src/physics/Diagnostics.hpp"
#include "../src/physics/UniformField.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

bool closeRelative(double actual, double expected, double tolerance = 1.0e-12) {
    return std::abs(actual - expected) <= tolerance * std::abs(expected);
}

void checkRelativity() {
    constexpr double energy_eV = 10000.0;
    const double energy_J = energy_eV * C::elementaryCharge_C;
    const double p = Relativity::momentumMagnitudeFromKineticEnergy_J(
        energy_J, C::electronMass_kg);
    const Vector3 momentum(p, 0.0, 0.0);
    const double gamma = Relativity::gammaFromMomentum(momentum, C::electronMass_kg);
    const double expectedGamma = 1.0 + energy_J /
        (C::electronMass_kg * C::speedOfLight_m_s * C::speedOfLight_m_s);

    require(closeRelative(gamma, expectedGamma), "incorrect Lorentz factor");
    require(closeRelative(Relativity::kineticEnergyFromMomentum_J(momentum, C::electronMass_kg),
                          energy_J), "energy and momentum do not round-trip");
    require(Relativity::velocityFromMomentum(momentum, C::electronMass_kg).norm() <
            C::speedOfLight_m_s, "particle velocity exceeds the speed of light");
    require(Relativity::kineticEnergyFromMomentum_J(Vector3{}, C::electronMass_kg) == 0.0,
            "resting particle has nonzero kinetic energy");
}

void checkDiagnostics() {
    constexpr double B_T = 5.0e-5;
    const double p = Relativity::momentumMagnitudeFromKineticEnergy_J(
        10000.0 * C::elementaryCharge_C, C::electronMass_kg);
    const Vector3 momentum(0.6 * p, 0.0, 0.8 * p);
    const Vector3 field(0.0, 0.0, B_T);
    const auto initial = evaluateDiagnostics(momentum, field, C::electronMass_kg);

    require(closeRelative(initial.energy_eV, 10000.0), "incorrect kinetic energy in eV");
    require(closeRelative(initial.pitchAngle_rad, std::acos(0.8)), "incorrect pitch angle");
    require(closeRelative(initial.parallelMomentum_kg_m_s, 0.8 * p),
            "incorrect parallel momentum");
    require(closeRelative(initial.magneticMoment_J_T,
                          (0.6 * p) * (0.6 * p) / (2.0 * C::electronMass_kg * B_T)),
            "incorrect magnetic moment");

    const UniformField uniformField(field);
    const BorisIntegrator integrator(uniformField, -C::elementaryCharge_C, C::electronMass_kg);
    const double period_s = 2.0 * std::numbers::pi *
        Relativity::gammaFromMomentum(momentum, C::electronMass_kg) * C::electronMass_kg /
        (C::elementaryCharge_C * B_T);
    constexpr int stepsPerPeriod = 100;
    ParticleState state{Vector3{}, momentum, 0.0};
    double maxEnergyError = 0.0;
    double maxPitchError = 0.0;
    double maxMomentError = 0.0;
    for (int i = 0; i < 2 * stepsPerPeriod; ++i) {
        state = integrator.step(state, period_s / stepsPerPeriod);
        const auto current = evaluateDiagnostics(
            state.p, uniformField.evaluate(state.r), C::electronMass_kg);
        maxEnergyError = std::max(maxEnergyError,
            std::abs(current.energy_eV / initial.energy_eV - 1.0));
        maxPitchError = std::max(maxPitchError,
            std::abs(current.pitchAngle_rad - initial.pitchAngle_rad));
        maxMomentError = std::max(maxMomentError,
            std::abs(current.magneticMoment_J_T / initial.magneticMoment_J_T - 1.0));
    }

    require(maxEnergyError < 1.0e-11, "energy drifted in uniform magnetic field");
    require(maxPitchError < 1.0e-11, "pitch angle drifted in uniform magnetic field");
    require(maxMomentError < 1.0e-11, "magnetic moment drifted in uniform magnetic field");

    bool zeroMomentumRejected = false;
    bool zeroFieldRejected = false;
    try { (void)evaluateDiagnostics(Vector3{}, field, C::electronMass_kg); }
    catch (const std::invalid_argument&) { zeroMomentumRejected = true; }
    try { (void)evaluateDiagnostics(momentum, Vector3{}, C::electronMass_kg); }
    catch (const std::invalid_argument&) { zeroFieldRejected = true; }
    require(zeroMomentumRejected && zeroFieldRejected,
            "undefined diagnostics were not rejected");

    std::cout << "Maximum relative energy error: " << maxEnergyError << '\n';
    std::cout << "Maximum pitch angle error (rad): " << maxPitchError << '\n';
    std::cout << "Maximum relative magnetic moment error: " << maxMomentError << '\n';
}

} // namespace

int main() {
    try {
        checkRelativity();
        checkDiagnostics();
    } catch (const std::exception& error) {
        std::cerr << "Diagnostics check failed: " << error.what() << '\n';
        return 1;
    }
}
