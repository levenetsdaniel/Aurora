#include "../src/simulation/SimulationEngine.hpp"
#include "../src/physics/UniformField.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

template <typename Action>
void requireInvalidArgument(Action action, const char* message) {
    try {
        action();
    } catch (const std::invalid_argument&) {
        return;
    }
    throw std::runtime_error(message);
}

constexpr double startRadius_m = C::earthRadius_m + 200000.0;
constexpr double atmosphereRadius_m = C::earthRadius_m + 150000.0;
constexpr double escapeRadius_m = C::earthRadius_m + 300000.0;
constexpr double step_s = 1.0e-5;

SimulationConfig baseConfig() {
    return {step_s, 150000.0, escapeRadius_m, 0.01, 10000};
}

double electronMomentum() {
    return Relativity::momentumMagnitudeFromKineticEnergy_J(
        1000.0 * C::elementaryCharge_C, C::electronMass_kg);
}

ParticleState radialState(double direction) {
    return {Vector3(startRadius_m, 0.0, 0.0),
            Vector3(direction * electronMomentum(), 0.0, 0.0), 0.0};
}

void checkVectorFinite() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    require(Vector3(1.0, -2.0, 3.0).isFinite(), "finite vector rejected");
    require(!Vector3(nan, 0.0, 0.0).isFinite(), "NaN vector accepted");
    require(!Vector3(0.0, inf, 0.0).isFinite(), "infinite vector accepted");
}

void checkBoundaryOutcomes(const UniformField& field) {
    const SimulationEngine engine(field, -C::elementaryCharge_C,
                                  C::electronMass_kg, baseConfig());
    const double speed_m_s = Relativity::velocityFromMomentum(
        radialState(1.0).p, C::electronMass_kg).norm();

    const auto inward = engine.run(radialState(-1.0));
    const double entryTime_s = (startRadius_m - atmosphereRadius_m) / speed_m_s;
    require(inward.outcome == ParticleOutcome::Precipitated, "inward particle did not precipitate");
    require(inward.steps > 0 && inward.steps < baseConfig().maxSteps, "invalid precipitation step count");
    require(std::abs(inward.eventState.r.norm() - atmosphereRadius_m) < 1e-3,
            "precipitation state was not localized to the atmospheric boundary");
    require(inward.eventState.time_s >= entryTime_s - 1e-12 &&
            inward.eventState.time_s <= entryTime_s + step_s + 1e-12,
            "precipitation time disagrees with straight-line motion");

    const auto outward = engine.run(radialState(1.0));
    const double escapeTime_s = (escapeRadius_m - startRadius_m) / speed_m_s;
    require(outward.outcome == ParticleOutcome::Escaped, "outward particle did not escape");
    require(std::abs(outward.eventState.r.norm() - escapeRadius_m) < 1e-3,
            "escape state was not localized to the outer boundary");
    require(outward.eventState.time_s >= escapeTime_s - 1e-12 &&
            outward.eventState.time_s <= escapeTime_s + step_s + 1e-12,
            "escape time disagrees with straight-line motion");

    const auto again = engine.run(radialState(-1.0));
    require(again.outcome == inward.outcome && again.steps == inward.steps &&
            again.eventState.time_s == inward.eventState.time_s &&
            again.eventState.r.x == inward.eventState.r.x,
            "identical runs gave different results");
}

void checkLimits(const UniformField& field) {
    SimulationConfig timeLimited = baseConfig();
    timeLimited.maxTime_s = 0.001003; // Deliberately not a multiple of step_s.
    const SimulationEngine timeEngine(field, -C::elementaryCharge_C,
                                      C::electronMass_kg, timeLimited);
    const auto byTime = timeEngine.run(radialState(-1.0));
    require(byTime.outcome == ParticleOutcome::Unresolved, "time limit did not stop the run");
    require(std::abs(byTime.eventState.time_s - timeLimited.maxTime_s) < 1e-12,
            "time-limited result has the wrong final time");

    SimulationConfig stepLimited = baseConfig();
    stepLimited.maxSteps = 5;
    const SimulationEngine stepEngine(field, -C::elementaryCharge_C,
                                      C::electronMass_kg, stepLimited);
    const auto bySteps = stepEngine.run(radialState(-1.0));
    require(bySteps.outcome == ParticleOutcome::Unresolved, "step limit did not stop the run");
    require(bySteps.steps == 5 && std::abs(bySteps.eventState.time_s - 5 * step_s) < 1e-12,
            "step-limited result has incorrect steps or time");
}

void checkInvalidInputs(const UniformField& field) {
    SimulationConfig config = baseConfig();
    config.dt_s = 0.0;
    requireInvalidArgument([&] {
        SimulationEngine engine(field, -C::elementaryCharge_C, C::electronMass_kg, config);
    }, "zero time step accepted");

    config = baseConfig();
    config.maxSteps = 0;
    requireInvalidArgument([&] {
        SimulationEngine engine(field, -C::elementaryCharge_C, C::electronMass_kg, config);
    }, "zero step limit accepted");

    config = baseConfig();
    config.escapeRadius_m = atmosphereRadius_m;
    requireInvalidArgument([&] {
        SimulationEngine engine(field, -C::elementaryCharge_C, C::electronMass_kg, config);
    }, "outer boundary inside atmosphere accepted");

    const SimulationEngine engine(field, -C::elementaryCharge_C,
                                  C::electronMass_kg, baseConfig());
    ParticleState state = radialState(-1.0);
    state.r.x = std::numeric_limits<double>::quiet_NaN();
    requireInvalidArgument([&] { (void)engine.run(state); }, "NaN position accepted");

    state = radialState(-1.0);
    state.p.y = std::numeric_limits<double>::infinity();
    requireInvalidArgument([&] { (void)engine.run(state); }, "infinite momentum accepted");

    state = radialState(-1.0);
    state.time_s = -1.0;
    requireInvalidArgument([&] { (void)engine.run(state); }, "negative initial time accepted");

    state = radialState(-1.0);
    state.r.x = atmosphereRadius_m - 1.0;
    requireInvalidArgument([&] { (void)engine.run(state); }, "start inside atmosphere accepted");

    state = radialState(1.0);
    state.r.x = escapeRadius_m + 1.0;
    requireInvalidArgument([&] { (void)engine.run(state); }, "start outside escape boundary accepted");

    state = radialState(-1.0);
    state.time_s = baseConfig().maxTime_s;
    requireInvalidArgument([&] { (void)engine.run(state); }, "initial time at limit accepted");
}

} // namespace

int main() {
    try {
        const UniformField field(Vector3(1e-9, 0.0, 0.0));
        checkVectorFinite();
        checkBoundaryOutcomes(field);
        checkLimits(field);
        checkInvalidInputs(field);
        std::cout << "SimulationEngine checks passed\n";
    } catch (const std::exception& error) {
        std::cerr << "SimulationEngine check failed: " << error.what() << '\n';
        return 1;
    }
}
