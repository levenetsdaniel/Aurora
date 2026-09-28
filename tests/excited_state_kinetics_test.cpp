#include "../src/emission/ExcitedStateKinetics.hpp"
#include "../src/io/KineticsCsvLoader.hpp"

#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>

#ifndef AURORA_SOURCE_DIR
#error AURORA_SOURCE_DIR is required for the kinetics data test
#endif

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

CrossSectionDatabase crossSections() {
    return CrossSectionDatabase({
        {"O-1D", "O", CollisionType::Excitation, "O(1D)", 2.0,
         TabulatedFunction({{2.0, 0.0}, {10.0, 1e-20}}), "synthetic test"}
    });
}

ExcitationModel production(double markerRate_s1) {
    const auto database = crossSections();
    ExcitationModel model(database,
        {{C::earthRadius_m + 100.0, -0.5, -0.5}, {1.0, 1.0, 1.0}, 1, 1, 1});
    if (markerRate_s1 > 0.0) {
        ElectronCascadeResult cascade;
        ElectronTrackRecord track;
        ElectronCollisionEvent event;
        event.position_m = {C::earthRadius_m + 100.5, 0.0, 0.0};
        event.channelId = "O-1D";
        event.type = CollisionType::Excitation;
        track.transport.events.push_back(event);
        cascade.tracks.push_back(track);
        model.accumulate(cascade, markerRate_s1);
    }
    return model;
}

KineticsCatalog syntheticCatalog() {
    return KineticsCatalog(
        {{"O(1D)", "O", 2.0, "synthetic test"}},
        {{"O(1D)", "O", 1e-17, 300.0, 0.0, 0.0, "synthetic test"},
         {"O(1D)", "O2", 0.0, 300.0, 0.0, 0.0, "synthetic zero"},
         {"O(1D)", "N2", 0.0, 300.0, 0.0, 0.0, "synthetic zero"}},
        {{{"O", "O(1D)"}, "O(1D)"}});
}

AtmosphereModel gas(double nO_m3) {
    return AtmosphereModel({
        {0.0, nO_m3, 0.0, 0.0, 200.0},
        {1000.0, nO_m3, 0.0, 0.0, 200.0}
    }, "synthetic atmosphere");
}

void checkPhysicalLimits() {
    const auto catalog = syntheticCatalog();
    const auto excited = production(8.0);
    const auto vacuum = gas(0.0);
    const ExcitedStateKinetics vacuumKinetics(catalog, vacuum);
    const auto noQuenching = vacuumKinetics.solve(excited, {"O", "O(1D)"}, 0, 0, 0);
    require(near(noQuenching.production_m3_s1, 8.0) &&
            near(noQuenching.radiativeDecay_s1, 2.0) &&
            near(noQuenching.quenching_s1, 0.0) &&
            near(noQuenching.effectiveLifetime_s, 0.5) &&
            near(noQuenching.density_m3, 4.0),
            "zero-quenching limit P/A is wrong");

    const auto denseGas = gas(1e18);
    const ExcitedStateKinetics denseKinetics(catalog, denseGas);
    const auto quenched = denseKinetics.solve(excited, {"O", "O(1D)"}, 0, 0, 0);
    require(near(quenched.quenching_s1, 10.0) &&
            near(quenched.effectiveLifetime_s, 1.0 / 12.0) &&
            near(quenched.density_m3, 8.0 / 12.0) &&
            quenched.density_m3 < noQuenching.density_m3,
            "strong quenching did not suppress the population");

    const auto empty = production(0.0);
    const auto noProduction = denseKinetics.solve(empty, {"O", "O(1D)"}, 0, 0, 0);
    require(noProduction.production_m3_s1 == 0.0 &&
            noProduction.density_m3 == 0.0,
            "zero production generated an excited population");
    require(denseKinetics.solveGrid(excited, {"O", "O(1D)"}).size() == 1,
            "grid solution has wrong number of cells");
}

void checkPublishedTables() {
    const auto catalog = KineticsCsvLoader::load(
        std::string(AURORA_SOURCE_DIR) + "/data/kinetics");
    const auto& red = catalog.resolve({"O", "O(1D)"});
    const auto& green = catalog.resolve({"O", "O -> O*((2p4)1S)"});
    require(near(red.totalDecay_s1, 7.45e-3) &&
            near(green.totalDecay_s1, 1.3354) &&
            near(catalog.radiativeLifetime_s({"O", "O(1D)"}), 1.0 / 7.45e-3) &&
            !red.source.empty() && !green.source.empty() &&
            catalog.termsFor({"O", "O(1D)"}).size() == 4,
            "published state rates, aliases or provenance were not loaded");
    const AtmosphereSample sample{150.0, 1e18, 1e18, 1e18, 200.0};
    const double expectedRed = 9.2e-19 * 1e18 +
        (1e-18 + 2.6e-17 * std::exp(67.0 / 200.0)) * 1e18 +
        2.3e-17 * 1e18;
    require(near(catalog.quenchingFrequency_s1("O(1D)", sample), expectedRed),
            "published O(1D) quenching terms or cm3-to-m3 conversion are wrong");
    const double expectedGreen = 5e-18 * std::exp(-301.0 / 200.0) * 1e18 +
        4e-18 * std::exp(-865.0 / 200.0) * 1e18 + 1e-23 * 1e18;
    require(near(catalog.quenchingFrequency_s1("O(1S)", sample), expectedGreen),
            "published O(1S) temperature dependence is wrong");
    const auto excited = production(8.0);
    const auto atmosphere = gas(1e18);
    const ExcitedStateKinetics kinetics(catalog, atmosphere);
    const auto population = kinetics.solve(excited, {"O", "O(1D)"}, 0, 0, 0);
    const AtmosphereSample localSample = atmosphere.sample(100.5, 0.0, 0.0, 0.0);
    require(near(population.density_m3,
                 8.0 / (red.totalDecay_s1 +
                        catalog.quenchingFrequency_s1(red.id, localSample))),
            "CSV-to-population integration is wrong");
    requireThrows<std::out_of_range>([&] {
        (void)catalog.resolve({"N2", "N2(A)"});
    }, "unknown state was silently given oxygen kinetics");
}

void checkValidation() {
    const auto database = crossSections();
    ExcitationModel excitation(database,
        {{C::earthRadius_m + 100.0, -0.5, -0.5}, {1.0, 1.0, 1.0}, 1, 1, 1});
    const auto catalog = syntheticCatalog();
    const auto atmosphere = gas(0.0);
    const ExcitedStateKinetics kinetics(catalog, atmosphere);
    requireThrows<std::out_of_range>([&] {
        (void)kinetics.solve(excitation, {"O", "missing"}, 0, 0, 0);
    }, "missing kinetic state was accepted");
    requireThrows<std::invalid_argument>([] {
        (void)KineticsCatalog(
            {{"bad", "O", 0.0, "source"}},
            {}, {{{"O", "bad"}, "bad"}});
    }, "zero radiative decay rate was accepted");
    requireThrows<std::invalid_argument>([] {
        (void)KineticsCatalog(
            {{"test", "O", 1.0, "source"}},
            {{"test", "O", 1e-17, 300.0, 0.0, 0.0, "source"}},
            {{{"O", "test"}, "test"}});
    }, "missing O2 and N2 quenching data were silently treated as zero");
    requireThrows<std::invalid_argument>([&] {
        (void)catalog.quenchingFrequency_s1("O(1D)",
            {100.0, -1.0, 0.0, 0.0, 200.0});
    }, "negative atmosphere density was accepted");
    std::istringstream badStates("id,species,A_total_s1,source\nO(1D),O,-1,source\n");
    std::istringstream terms("id,quencher,k0_m3_s,T_ref_K,beta,activation_K,source\n"
                             "O(1D),O,0,300,0,0,source\n"
                             "O(1D),O2,0,300,0,0,source\n"
                             "O(1D),N2,0,300,0,0,source\n");
    std::istringstream aliases("species,event_state,id\nO,O(1D),O(1D)\n");
    requireThrows<std::invalid_argument>([&] {
        (void)KineticsCsvLoader::load(badStates, terms, aliases);
    }, "negative coefficient in CSV was accepted");
}

} // namespace

int main() {
    try {
        checkPhysicalLimits();
        checkPublishedTables();
        checkValidation();
        std::cout << "ExcitedStateKinetics checks passed\n";
    } catch (const std::exception& error) {
        std::cerr << "ExcitedStateKinetics check failed: " << error.what() << '\n';
        return 1;
    }
}
