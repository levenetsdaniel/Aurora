#include "../src/numerics/BorisIntegrator.hpp"
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

struct OrbitResult {
    double positionError_m;
    double relativeMomentumError;
};

OrbitResult runUniformOrbit(int stepsPerPeriod) {
    constexpr double field_T = 5.0e-5;
    constexpr double kineticEnergy_eV = 10000.0;
    constexpr double perpendicularFraction = 0.8;
    constexpr double periods = 3.25;

    const UniformField field(Vector3(0.0, 0.0, field_T));
    const double charge_C = -C::elementaryCharge_C;
    const double mass_kg = C::electronMass_kg;
    const double gamma = 1.0 + kineticEnergy_eV * C::elementaryCharge_C /
        (mass_kg * C::speedOfLight_m_s * C::speedOfLight_m_s);
    const double momentum = mass_kg * C::speedOfLight_m_s * std::sqrt(gamma * gamma - 1.0);
    const double pPerp = momentum * perpendicularFraction;
    const double pParallel = momentum * std::sqrt(1.0 - perpendicularFraction * perpendicularFraction);
    const double omega = charge_C * field_T / (gamma * mass_kg);
    const double period_s = 2.0 * std::numbers::pi / std::abs(omega);
    const double dt_s = period_s / stepsPerPeriod;
    const int steps = static_cast<int>(periods * stepsPerPeriod);

    BorisIntegrator integrator(field, charge_C, mass_kg);
    ParticleState state{Vector3{}, Vector3(pPerp, 0.0, pParallel), 0.0};
    double maxRelativeMomentumError = 0.0;
    for (int i = 0; i < steps; ++i) {
        state = integrator.step(state, dt_s);
        maxRelativeMomentumError = std::max(maxRelativeMomentumError,
            std::abs(state.p.norm() - momentum) / momentum);
    }

    const double time_s = steps * dt_s;
    const double vPerp = pPerp / (gamma * mass_kg);
    const double vParallel = pParallel / (gamma * mass_kg);
    const Vector3 exactPosition(
        vPerp * std::sin(omega * time_s) / omega,
        vPerp * (std::cos(omega * time_s) - 1.0) / omega,
        vParallel * time_s);

    require(std::abs(state.time_s - time_s) < 1.0e-12 * period_s,
            "incorrect elapsed time");
    return {(state.r - exactPosition).norm(), maxRelativeMomentumError};
}

} // namespace

int main() {
    try {
        const auto coarse = runUniformOrbit(50);
        const auto medium = runUniformOrbit(100);
        const auto fine = runUniformOrbit(200);
        const auto finest = runUniformOrbit(400);

        require(medium.positionError_m < coarse.positionError_m / 3.0,
                "orbit error did not converge from 50 to 100 steps per period");
        require(fine.positionError_m < medium.positionError_m / 3.0,
                "orbit error did not converge from 100 to 200 steps per period");
        require(finest.positionError_m < fine.positionError_m / 3.0,
                "orbit error did not converge from 200 to 400 steps per period");
        require(finest.relativeMomentumError < 1.0e-11,
                "momentum magnitude drifted in a magnetic-only orbit");

        std::cout << "Position errors (m), 50/100/200/400 steps per period: "
                  << coarse.positionError_m << ", " << medium.positionError_m << ", "
                  << fine.positionError_m << ", " << finest.positionError_m << '\n';
        std::cout << "Maximum relative momentum error: "
                  << finest.relativeMomentumError << '\n';
    } catch (const std::exception& error) {
        std::cerr << "Boris integrator check failed: " << error.what() << '\n';
        return 1;
    }
}
