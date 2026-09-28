#include "../src/atmosphere/ElectronCascade.hpp"

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

bool near(double actual, double expected) {
    return std::abs(actual - expected) < 1e-10;
}

AtmosphereModel gas(double density_m3) {
    return AtmosphereModel({
        {0.0, density_m3, 0.0, 0.0, 200.0},
        {2000.0, density_m3, 0.0, 0.0, 200.0}
    }, "synthetic homogeneous oxygen gas");
}

CrossSectionDatabase excitation() {
    return CrossSectionDatabase({
        {"O-excitation", "O", CollisionType::Excitation, "test state", 5.0,
         TabulatedFunction({{5.0, 0.0}, {6.0, 1e-20}, {20.0, 1e-20}}),
         "synthetic test"}
    });
}

CrossSectionDatabase ionization() {
    return CrossSectionDatabase({
        {"O-ionization", "O", CollisionType::Ionization, "", 5.0,
         TabulatedFunction({{5.0, 0.0}, {6.0, 1e-20}, {20.0, 1e-20}}),
         "synthetic test"}
    });
}

ElectronCascadeConfig config() {
    return {100.0, 5.0, 100, 100, 1.0};
}

void checkNoCollisionAndEscape() {
    const auto emptyGas = gas(0.0);
    const auto crossSections = excitation();
    const ElectronTransportMC transport(emptyGas, crossSections);
    const ElectronCascade cascade(transport, config());
    Random rng(11);
    const auto start = Vector3(C::earthRadius_m + 1000.0, 0.0, 0.0);
    const auto noCollision = cascade.run(start, {1.0, 0.0, 0.0}, 10.0, rng);
    require(noCollision.tracks.size() == 1 &&
            noCollision.tracks[0].transport.events.empty() &&
            near(noCollision.ledger.remaining_eV, 10.0) &&
            near(noCollision.ledger.residual_eV(), 0.0) &&
            !noCollision.limitReached,
            "collision-free track did not retain its energy");

    const auto escaped = cascade.run({C::earthRadius_m + 1990.0, 0.0, 0.0},
                                     {1.0, 0.0, 0.0}, 10.0, rng);
    require(escaped.tracks.size() == 1 &&
            escaped.tracks[0].transport.exitedAtmosphere &&
            near(escaped.ledger.escaped_eV, 10.0) &&
            near(escaped.ledger.remaining_eV, 0.0) &&
            near(escaped.ledger.residual_eV(), 0.0),
            "escaping track was not recorded as outgoing energy");
}

void checkOneExcitationAndThermalization() {
    const auto denseGas = gas(1e20);
    const auto crossSections = excitation();
    const ElectronTransportMC transport(denseGas, crossSections);
    const ElectronCascade cascade(transport, config());
    Random rng(12);
    const auto start = Vector3(C::earthRadius_m + 1000.0, 0.0, 0.0);
    const auto result = cascade.run(start, {1.0, 0.0, 0.0}, 10.0, rng);
    require(result.tracks.size() == 1 && result.tracks[0].transport.events.size() == 1 &&
            near(result.ledger.excitation_eV, 5.0) &&
            near(result.ledger.remaining_eV, 5.0) &&
            near(result.ledger.residual_eV(), 0.0) &&
            !result.limitReached,
            "single excitation was not accounted for");

    const auto cold = cascade.run(start, {1.0, 0.0, 0.0}, 0.5, rng);
    require(cold.tracks.size() == 1 && cold.tracks[0].thermalizedWithoutTransport &&
            near(cold.ledger.thermal_eV, 0.5) &&
            near(cold.ledger.residual_eV(), 0.0),
            "sub-thermal primary electron was not deposited");
}

void checkSecondaryQueueAndLimits() {
    const auto denseGas = gas(1e20);
    const auto crossSections = ionization();
    const ElectronTransportMC transport(denseGas, crossSections);
    const ElectronCascade cascade(transport, config());
    Random rng(13);
    const auto start = Vector3(C::earthRadius_m + 1000.0, 0.0, 0.0);
    const auto result = cascade.run(start, {1.0, 0.0, 0.0}, 10.0, rng);
    require(result.tracks.size() == 2 &&
            result.tracks[1].parentIndex == 0 && result.tracks[1].generation == 1 &&
            result.tracks[0].ledger.secondaryTransfer_eV > 0.0 &&
            near(result.secondaryGenerated_eV,
                 result.tracks[1].initialEnergy_eV) &&
            near(result.ledger.ionization_eV, 5.0) &&
            near(result.ledger.thermal_eV + result.ledger.remaining_eV +
                 result.ledger.escaped_eV, 5.0) &&
            near(result.ledger.residual_eV(), 0.0) &&
            near(result.tracks[0].ledger.residual_eV(), 0.0) &&
            near(result.tracks[1].ledger.residual_eV(), 0.0) &&
            !result.limitReached,
            "secondary electron was not transported or energy was double counted");

    auto limitedConfig = config();
    limitedConfig.maxTracks = 1;
    const ElectronCascade limited(transport, limitedConfig);
    Random limitedRng(13);
    const auto truncated = limited.run(start, {1.0, 0.0, 0.0}, 10.0,
                                       limitedRng);
    require(truncated.tracks.size() == 1 && truncated.unprocessedTracks == 1 &&
            truncated.limitReached && near(truncated.ledger.residual_eV(), 0.0),
            "track cap lost the queued secondary's energy");
}

void checkInvalidInput() {
    const auto emptyGas = gas(0.0);
    const auto crossSections = excitation();
    const ElectronTransportMC transport(emptyGas, crossSections);
    auto badConfig = config();
    badConfig.maxTracks = 0;
    requireThrows<std::invalid_argument>([&] {
        (void)ElectronCascade(transport, badConfig);
    }, "zero track limit was accepted");
    const ElectronCascade cascade(transport, config());
    Random rng(14);
    requireThrows<std::invalid_argument>([&] {
        (void)cascade.run({C::earthRadius_m + 1000.0, 0.0, 0.0},
                          {}, 10.0, rng);
    }, "zero direction was accepted");
}

void checkEventLimitKeepsEnergy() {
    const auto denseGas = gas(1e20);
    const CrossSectionDatabase elastic({
        {"O-elastic", "O", CollisionType::Elastic, "", 0.0,
         TabulatedFunction({{0.0, 1e-20}, {20.0, 1e-20}}),
         "synthetic test"}
    });
    const ElectronTransportMC transport(denseGas, elastic);
    auto limitedConfig = config();
    limitedConfig.maxEventsPerTrack = 1;
    const ElectronCascade cascade(transport, limitedConfig);
    Random rng(15);
    const auto result = cascade.run({C::earthRadius_m + 1000.0, 0.0, 0.0},
                                    {1.0, 0.0, 0.0}, 10.0, rng);
    require(result.tracks.size() == 1 && result.tracks[0].transport.events.size() == 1 &&
            result.limitReached && near(result.ledger.remaining_eV, 10.0) &&
            near(result.ledger.residual_eV(), 0.0),
            "event cap lost the unfinished primary's energy");
}

} // namespace

int main() {
    try {
        checkNoCollisionAndEscape();
        checkOneExcitationAndThermalization();
        checkSecondaryQueueAndLimits();
        checkEventLimitKeepsEnergy();
        checkInvalidInput();
        std::cout << "ElectronCascade checks passed\n";
    } catch (const std::exception& error) {
        std::cerr << "ElectronCascade check failed: " << error.what() << '\n';
        return 1;
    }
}
