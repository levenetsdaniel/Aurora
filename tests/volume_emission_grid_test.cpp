#include "../src/emission/VolumeEmissionGrid.hpp"
#include "../src/io/KineticsCsvLoader.hpp"
#include "../src/io/SpectralLineCsvLoader.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>

#ifndef AURORA_SOURCE_DIR
#error AURORA_SOURCE_DIR is required
#endif

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

template <typename Exception, typename Action>
void requireThrows(Action action, const char* message) {
    try { action(); } catch (const Exception&) { return; }
    throw std::runtime_error(message);
}

bool near(double actual, double expected) {
    return std::abs(actual - expected) <= 1e-12 * std::max(1.0, std::abs(expected));
}

KineticsCatalog kineticsCatalog() {
    return KineticsCatalog(
        {{"red", "O", 2.0, "test"}, {"green", "O", 4.0, "test"}},
        {{"red", "O", 1e-18, 300, 0, 0, "test"},
         {"red", "O2", 0, 300, 0, 0, "test"},
         {"red", "N2", 0, 300, 0, 0, "test"},
         {"green", "O", 0, 300, 0, 0, "test"},
         {"green", "O2", 0, 300, 0, 0, "test"},
         {"green", "N2", 0, 300, 0, 0, "test"}},
        {{{"O", "red-event"}, "red"},
         {{"O", "red-alias"}, "red"},
         {{"O", "green-event"}, "green"}});
}

CrossSectionDatabase crossSections() {
    return CrossSectionDatabase({
        {"red", "O", CollisionType::Excitation, "red-event", 1,
         TabulatedFunction({{1, 0}, {10, 1e-20}}), "test"},
        {"red-alias", "O", CollisionType::Excitation, "red-alias", 1,
         TabulatedFunction({{1, 0}, {10, 1e-20}}), "test"},
        {"green", "O", CollisionType::Excitation, "green-event", 1,
         TabulatedFunction({{1, 0}, {10, 1e-20}}), "test"}
    });
}

ExcitationModel excitation(const CrossSectionDatabase& database, bool filled) {
    ExcitationModel model(database,
        {{C::earthRadius_m + 100, -0.5, -0.5}, {1, 1, 1}, 1, 1, 1});
    if (!filled) return model;
    ElectronCascadeResult cascade;
    ElectronTrackRecord track;
    for (const char* id : {"red", "red-alias", "green"}) {
        ElectronCollisionEvent event;
        event.channelId = id;
        event.position_m = {C::earthRadius_m + 100.5, 0, 0};
        event.type = CollisionType::Excitation;
        track.transport.events.push_back(event);
    }
    cascade.tracks.push_back(track);
    model.accumulate(cascade, 8.0);
    return model;
}

void checkEmission() {
    const auto catalog = kineticsCatalog();
    const SpectralLineCatalog lines(
        {{"O_6300", 630.0, "red", 1.0, "test"},
         {"O_5577", 557.7, "green", 3.0, "test"}}, catalog);
    const AtmosphereModel atmosphere({{0, 1e18, 0, 0, 300},
                                      {1000, 1e18, 0, 0, 300}}, "test");
    const ExcitedStateKinetics kinetics(catalog, atmosphere);
    const auto database = crossSections();
    const auto produced = excitation(database, true);
    const auto emission = VolumeEmissionGrid::fromKinetics(
        produced, kinetics, catalog, lines);
    require(near(emission.emissivity("O_6300", 0, 0, 0), 16.0 / 3.0),
            "red line did not sum aliases or apply quenching");
    require(near(emission.emissivity("O_5577", 0, 0, 0), 6.0),
            "green line did not use its own A_lambda");
    const auto empty = VolumeEmissionGrid::fromKinetics(
        excitation(database, false), kinetics, catalog, lines);
    require(empty.emissivity("O_6300", 0, 0, 0) == 0.0 &&
            empty.emissivity("O_5577", 0, 0, 0) == 0.0,
            "zero production emitted photons");
}

void checkPublishedLines() {
    const auto catalog = KineticsCsvLoader::load(
        std::string(AURORA_SOURCE_DIR) + "/data/kinetics");
    const auto lines = SpectralLineCsvLoader::load(
        std::string(AURORA_SOURCE_DIR) + "/data/kinetics/lines.csv", catalog);
    require(lines.lines().size() == 2 &&
            near(lines.lines()[0].transitionRate_s1, 1.26) &&
            near(lines.lines()[1].transitionRate_s1, 5.63e-3),
            "published line coefficients were not loaded");
    requireThrows<std::invalid_argument>([&] {
        (void)SpectralLineCatalog(
            {{"bad", 630, "O(1D)", 7.45e-3, "test"},
             {"bad2", 636, "O(1D)", 1.82e-3, "test"}}, catalog);
    }, "line coefficients exceeding A_total were accepted");
    std::istringstream malformed(
        "id,wavelength_nm,upper_state_id,A_line_s1,source\n"
        "bad,630,O(1D),nan,test\n");
    requireThrows<std::invalid_argument>([&] {
        (void)SpectralLineCsvLoader::load(malformed, catalog);
    }, "nonfinite line coefficient was accepted");
}

} // namespace

int main() {
    try {
        checkEmission();
        checkPublishedLines();
        std::cout << "VolumeEmissionGrid checks passed\n";
    } catch (const std::exception& error) {
        std::cerr << "VolumeEmissionGrid check failed: " << error.what() << '\n';
        return 1;
    }
}
