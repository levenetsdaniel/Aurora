#include "../src/simulation/LossCone.hpp"
#include "../src/physics/Constants.hpp"
#include "../src/physics/DipoleField.hpp"
#include "../src/physics/Relativity.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

template <typename Exception, typename Action>
void requireThrows(Action action, const char* message) {
    try {
        action();
    } catch (const Exception&) {
        return;
    }
    throw std::runtime_error(message);
}

// A deterministic test field: the component along -y grows on approach to Earth.
// With zero charge, pitch alone determines where p dot B changes sign.
class TurningField final : public MagneticField {
private:
    double startX_m;
    double scale_m;

public:
    TurningField(double startX_m, double scale_m)
        : startX_m(startX_m), scale_m(scale_m) {}

    Vector3 evaluate(const Vector3& r) const override {
        return Vector3(-1.0, -(startX_m - r.x) / scale_m, 0.0) * 1e-9;
    }
};

constexpr double startX_m = C::earthRadius_m + 200000.0;

ParticleState templateState() {
    const double momentum = Relativity::momentumMagnitudeFromKineticEnergy_J(
        1000.0 * C::elementaryCharge_C, C::electronMass_kg);
    return {Vector3(startX_m, 0.0, 0.0),
            Vector3(-std::cos(0.4), std::sin(0.4), 0.0) * momentum, 0.0};
}

SimulationConfig config() {
    return {1e-5, 150000.0, C::earthRadius_m + 300000.0, 0.02, 3000};
}

void checkRefinement() {
    const TurningField field(startX_m, 20000.0);
    const SimulationEngine engine(field, 0.0, C::electronMass_kg, config());
    const LossConeScanner scanner(engine);
    const ParticleState source = templateState();

    const LossConeBracket bracket = scanner.refine(source, 0.1, 0.8, 1e-4);
    require(bracket.precipitatingPitch_rad < bracket.mirroringPitch_rad &&
            bracket.mirroringPitch_rad - bracket.precipitatingPitch_rad <= 1e-4,
            "loss cone bracket was not refined to tolerance");
    require(bracket.iterations > 0, "scanner stopped without refining");
    require(bracket.precipitatingPitch_rad > 0.25 && bracket.mirroringPitch_rad < 0.5,
            "loss cone boundary is outside the expected range");
    const double momentum = source.p.norm();
    const Vector3 along(-1.0, 0.0, 0.0);
    const Vector3 perpendicular(0.0, 1.0, 0.0);
    const auto outcomeAt = [&](double pitch_rad) {
        ParticleState trial = source;
        trial.p = (along * std::cos(pitch_rad) +
                   perpendicular * std::sin(pitch_rad)) * momentum;
        require(std::abs(trial.p.norm() / momentum - 1.0) < 1e-12,
                "changing pitch changed momentum magnitude");
        return engine.run(trial).outcome;
    };
    require(outcomeAt(bracket.precipitatingPitch_rad) == ParticleOutcome::Precipitated &&
            outcomeAt(bracket.mirroringPitch_rad) == ParticleOutcome::Mirrored,
            "returned endpoints do not retain their outcomes");
    require(source.p.x == templateState().p.x && source.p.y == templateState().p.y,
            "scanner changed the template");

    const LossConeBracket repeated = scanner.refine(source, 0.1, 0.8, 1e-4);
    require(repeated.precipitatingPitch_rad == bracket.precipitatingPitch_rad &&
            repeated.mirroringPitch_rad == bracket.mirroringPitch_rad &&
            repeated.iterations == bracket.iterations,
            "identical scans produced different results");

    const LossConeBracket alreadyNarrow = scanner.refine(source, 0.1, 0.8, 0.71);
    require(alreadyNarrow.iterations == 0 &&
            alreadyNarrow.precipitatingPitch_rad == 0.1 &&
            alreadyNarrow.mirroringPitch_rad == 0.8,
            "scanner did not preserve an already narrow bracket");

    requireThrows<std::invalid_argument>([&] {
        (void)scanner.refine(source, 0.8, 0.1, 1e-4);
    }, "reversed pitch bracket was accepted");
    requireThrows<std::invalid_argument>([&] {
        (void)scanner.refine(source, 0.1, 0.2, 1e-4);
    }, "non-mirroring upper endpoint was accepted");
    requireThrows<std::runtime_error>([&] {
        (void)LossConeScanner(engine, 1).refine(source, 0.1, 0.8, 1e-6);
    }, "iteration limit was ignored");

    SimulationConfig shortConfig = config();
    shortConfig.maxTime_s = 1e-5;
    const SimulationEngine shortEngine(field, 0.0, C::electronMass_kg, shortConfig);
    requireThrows<std::invalid_argument>([&] {
        (void)LossConeScanner(shortEngine).refine(source, 0.1, 0.8, 1e-4);
    }, "unresolved endpoint was treated as a valid bracket");

    SimulationConfig midpointConfig = config();
    midpointConfig.maxTime_s = 0.00275;
    const SimulationEngine midpointEngine(field, 0.0, C::electronMass_kg, midpointConfig);
    requireThrows<std::runtime_error>([&] {
        (void)LossConeScanner(midpointEngine).refine(source, 0.1, 0.8, 1e-4);
    }, "unresolved midpoint was silently classified");
}

void checkDipoleAgainstAdiabaticEstimate() {
    const DipoleField field(1e-7);
    const double equatorialRadius_m = 2.0 * C::earthRadius_m;
    const double atmosphereRadius_m = C::earthRadius_m + 150000.0;
    const double cosineSquared = atmosphereRadius_m / equatorialRadius_m;
    const double sine = std::sqrt(1.0 - cosineSquared);
    const double cosine = std::sqrt(cosineSquared);
    const Vector3 atmospherePoint(atmosphereRadius_m * cosine, 0.0,
                                  atmosphereRadius_m * sine);
    const double equatorialB = field.evaluate(Vector3(equatorialRadius_m, 0.0, 0.0)).norm();
    const double atmosphereB = field.evaluate(atmospherePoint).norm();
    const double adiabaticPitch_rad = std::asin(std::sqrt(equatorialB / atmosphereB));

    const SimulationConfig physicalConfig{2e-5, 150000.0,
        3.0 * C::earthRadius_m, 2.0, 100000};
    const SimulationEngine engine(field, -C::elementaryCharge_C,
                                  C::electronMass_kg, physicalConfig);
    const double momentum = templateState().p.norm();
    const ParticleState source{Vector3(equatorialRadius_m, 0.0, 0.0),
        Vector3(std::sin(0.4), 0.0, std::cos(0.4)) * momentum, 0.0};

    const LossConeBracket bracket = LossConeScanner(engine).refine(source, 0.1, 0.6, 0.01);
    const double numericalPitch_rad =
        (bracket.precipitatingPitch_rad + bracket.mirroringPitch_rad) * 0.5;
    require(std::abs(numericalPitch_rad - adiabaticPitch_rad) < 0.1,
            "dipole loss cone differs substantially from the adiabatic estimate");

    SimulationConfig finerConfig = physicalConfig;
    finerConfig.dt_s *= 0.5;
    finerConfig.maxSteps *= 2;
    const SimulationEngine finerEngine(field, -C::elementaryCharge_C,
                                       C::electronMass_kg, finerConfig);
    const LossConeBracket finerBracket = LossConeScanner(finerEngine).refine(
        source, 0.1, 0.6, 0.01);
    const double finerPitch_rad =
        (finerBracket.precipitatingPitch_rad + finerBracket.mirroringPitch_rad) * 0.5;
    require(std::abs(finerPitch_rad - numericalPitch_rad) < 0.05,
            "dipole loss cone changed too much when time step was halved");

    ParticleState rotatedPhase = source;
    rotatedPhase.p = Vector3(0.0, std::sin(0.4), std::cos(0.4)) * momentum;
    const LossConeBracket phaseBracket = LossConeScanner(engine).refine(
        rotatedPhase, 0.1, 0.6, 0.01);
    const double phasePitch_rad =
        (phaseBracket.precipitatingPitch_rad + phaseBracket.mirroringPitch_rad) * 0.5;
    require(std::abs(phasePitch_rad - numericalPitch_rad) < 0.1,
            "dipole loss cone changed too much with gyro phase");
}

} // namespace

int main() {
    try {
        checkRefinement();
        checkDipoleAgainstAdiabaticEstimate();
        std::cout << "LossConeScanner checks passed\n";
    } catch (const std::exception& error) {
        std::cerr << "LossConeScanner check failed: " << error.what() << '\n';
        return 1;
    }
}
