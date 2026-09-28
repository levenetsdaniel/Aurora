#include "../src/emission/ExcitationModel.hpp"

#include <cmath>
#include <iostream>
#include <limits>
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
    return std::abs(actual - expected) <= 1e-12 * std::max(1.0, std::abs(expected));
}

CrossSectionDatabase channels() {
    return CrossSectionDatabase({
        {"O-1D", "O", CollisionType::Excitation, "1D", 2.0,
         TabulatedFunction({{2.0, 0.0}, {10.0, 1e-20}}), "synthetic test"},
        {"O2-1D", "O2", CollisionType::Excitation, "1D", 2.0,
         TabulatedFunction({{2.0, 0.0}, {10.0, 1e-20}}), "synthetic test"},
        {"O-ion", "O", CollisionType::Ionization, "", 5.0,
         TabulatedFunction({{5.0, 0.0}, {10.0, 1e-20}}), "synthetic test"}
    });
}

ElectronCollisionEvent event(Vector3 position_m, const std::string& id,
                             CollisionType type = CollisionType::Excitation) {
    ElectronCollisionEvent result;
    result.position_m = position_m;
    result.channelId = id;
    result.type = type;
    return result;
}

ElectronCascadeResult fixture() {
    ElectronCascadeResult result;
    ElectronTrackRecord primary;
    primary.transport.events.push_back(event({0.5, 0.5, 0.5}, "O-1D"));
    primary.transport.events.push_back(event({2.0, 0.5, 0.5}, "O2-1D"));
    primary.transport.events.push_back(event({0.5, 0.5, 0.5}, "O-ion",
                                             CollisionType::Ionization));
    ElectronTrackRecord secondary;
    secondary.parentIndex = 0;
    secondary.generation = 1;
    secondary.transport.events.push_back(event({0.6, 0.5, 0.5}, "O-1D"));
    result.tracks.push_back(std::move(primary));
    result.tracks.push_back(std::move(secondary));
    return result;
}

void checkAccumulationAndNormalization() {
    const auto crossSections = channels();
    const ExcitationGridConfig grid{{0.0, 0.0, 0.0}, {2.0, 2.0, 2.0}, 2, 1, 1};
    ExcitationModel model(crossSections, grid);
    const auto cascade = fixture();
    require(model.states().size() == 2,
            "states of different gas species were merged");
    require(model.accumulate(cascade, 4.0) == 3,
            "wrong number of excitation events accumulated");
    require(near(model.cellVolume(), 8.0) &&
            near(model.eventRate_s1({"O", "1D"}, 0, 0, 0), 8.0) &&
            near(model.productionRate_m3_s1({"O", "1D"}, 0, 0, 0), 1.0) &&
            near(model.productionRate_m3_s1({"O2", "1D"}, 1, 0, 0), 0.5) &&
            near(model.productionRate_m3_s1({"O", "1D"}, 1, 0, 0), 0.0),
            "species, state, cell or units were mixed up");

    ExcitationModel twoMarkers(crossSections, grid);
    twoMarkers.accumulate(cascade, 2.0);
    twoMarkers.accumulate(cascade, 2.0);
    require(near(twoMarkers.productionRate_m3_s1({"O", "1D"}, 0, 0, 0), 1.0),
            "splitting one physical rate between two markers changed production");
    twoMarkers.accumulate(cascade, 4.0);
    require(near(twoMarkers.productionRate_m3_s1({"O", "1D"}, 0, 0, 0), 2.0),
            "production did not scale linearly with marker rate");

    ExcitationModel largerCells(crossSections,
        {{0.0, 0.0, 0.0}, {4.0, 4.0, 4.0}, 1, 1, 1});
    largerCells.accumulate(cascade, 4.0);
    require(near(largerCells.productionRate_m3_s1({"O", "1D"}, 0, 0, 0), 0.125),
            "production rate did not scale inversely with cell volume");
}

void checkValidationAndOutsidePolicy() {
    const auto crossSections = channels();
    const ExcitationGridConfig grid{{0.0, 0.0, 0.0}, {2.0, 2.0, 2.0}, 2, 1, 1};
    ExcitationModel model(crossSections, grid);
    auto cascade = fixture();
    requireThrows<std::invalid_argument>([&] {
        model.accumulate(cascade, -1.0);
    }, "negative marker rate was accepted");
    requireThrows<std::invalid_argument>([&] {
        model.accumulate(cascade, std::numeric_limits<double>::infinity());
    }, "infinite marker rate was accepted");
    cascade.limitReached = true;
    requireThrows<std::invalid_argument>([&] {
        model.accumulate(cascade, 1.0);
    }, "incomplete cascade was accepted");
    cascade.limitReached = false;
    cascade.tracks[1].transport.events[0].position_m.x = 4.0;
    requireThrows<std::out_of_range>([&] {
        model.accumulate(cascade, 1.0);
    }, "upper half-open grid boundary was accepted");
    require(near(model.eventRate_s1({"O", "1D"}, 0, 0, 0), 0.0),
            "failed accumulation partially changed the grid");

    ExcitationModel ignoring(crossSections, grid, ExcitationOutOfBounds::Ignore);
    require(ignoring.accumulate(cascade, 3.0) == 2 &&
            near(ignoring.ignoredEventRate_s1(), 3.0),
            "ignored out-of-grid excitation was not reported");
    cascade.tracks[1].transport.events[0].channelId = "unknown";
    requireThrows<std::invalid_argument>([&] {
        ignoring.accumulate(cascade, 1.0);
    }, "unknown excitation channel was accepted");

    requireThrows<std::invalid_argument>([&] {
        (void)ExcitationModel(crossSections,
            {{0.0, 0.0, 0.0}, {0.0, 1.0, 1.0}, 1, 1, 1});
    }, "zero cell size was accepted");
    requireThrows<std::out_of_range>([&] {
        (void)model.productionRate_m3_s1({"O", "1D"}, 2, 0, 0);
    }, "out-of-grid cell index was accepted");
}

void checkRealCascadeInterface() {
    const CrossSectionDatabase crossSections({
        {"O-1D", "O", CollisionType::Excitation, "1D", 5.0,
         TabulatedFunction({{5.0, 0.0}, {6.0, 1e-20}, {20.0, 1e-20}}),
         "synthetic test"}
    });
    const AtmosphereModel gas({
        {0.0, 1e20, 0.0, 0.0, 200.0},
        {2000.0, 1e20, 0.0, 0.0, 200.0}
    }, "synthetic atmosphere");
    const ElectronTransportMC transport(gas, crossSections);
    const ElectronCascade cascade(transport, {100.0, 5.0, 10, 10, 1.0});
    Random rng(88);
    const auto result = cascade.run({C::earthRadius_m + 1000.0, 0.0, 0.0},
                                    {1.0, 0.0, 0.0}, 10.0, rng);
    require(!result.limitReached && result.tracks.size() == 1 &&
            result.tracks[0].transport.events.size() == 1,
            "cascade fixture did not generate exactly one excitation");
    ExcitationModel model(crossSections,
        {{C::earthRadius_m + 900.0, -100.0, -100.0},
         {200.0, 200.0, 200.0}, 1, 1, 1});
    require(model.accumulate(result, 5.0) == 1 &&
            near(model.productionRate_m3_s1({"O", "1D"}, 0, 0, 0),
                 5.0 / (200.0 * 200.0 * 200.0)),
            "real cascade event did not produce the expected cell rate");
}

} // namespace

int main() {
    try {
        checkAccumulationAndNormalization();
        checkValidationAndOutsidePolicy();
        checkRealCascadeInterface();
        std::cout << "ExcitationModel checks passed\n";
    } catch (const std::exception& error) {
        std::cerr << "ExcitationModel check failed: " << error.what() << '\n';
        return 1;
    }
}
