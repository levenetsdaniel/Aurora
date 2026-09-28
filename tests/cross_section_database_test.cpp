#include "../src/atmosphere/CrossSectionDatabase.hpp"
#include "../src/atmosphere/ReferenceCrossSections.hpp"
#include "../src/io/LxCatCrossSectionLoader.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

bool near(double actual, double expected) {
    return std::abs(actual - expected) <= 1e-12 * expected;
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

std::vector<CollisionChannel> exampleChannels() {
    return {
        {"N2-elastic", "N2", CollisionType::Elastic, "", 0.0,
         TabulatedFunction({{0.0, 4e-20}, {10.0, 2e-20}, {20.0, 1e-20}}),
         "synthetic test table"},
        {"N2-excitation", "N2", CollisionType::Excitation, "N2(A)", 5.0,
         TabulatedFunction({{5.0, 0.0}, {10.0, 2e-21}, {20.0, 4e-21}}),
         "synthetic test table"},
        {"N2-ionization", "N2", CollisionType::Ionization, "N2+", 12.0,
         TabulatedFunction({{12.0, 0.0}, {20.0, 8e-21}}),
         "synthetic test table"},
        {"O-elastic", "O", CollisionType::Elastic, "", 0.0,
         TabulatedFunction({{0.0, 1e-20}, {20.0, 1e-20}}),
         "synthetic test table"}
    };
}

void checkInterpolationAndChannels() {
    const CrossSectionDatabase db(exampleChannels());
    const CollisionChannel& elastic = db.channelById("N2-elastic");
    const CollisionChannel& excitation = db.channelById("N2-excitation");
    const CollisionChannel& ionization = db.channelById("N2-ionization");

    require(elastic.source == "synthetic test table" &&
            elastic.sigma_m2.minEnergy_eV() == 0.0 &&
            elastic.sigma_m2.maxEnergy_eV() == 20.0,
            "cross-section provenance or energy range was lost");
    require(near(db.sigma(elastic, 0.0), 4e-20) &&
            near(db.sigma(elastic, 5.0), 3e-20) &&
            near(db.sigma(elastic, 10.0), 2e-20) &&
            near(db.sigma(elastic, 20.0), 1e-20),
            "elastic cross-section interpolation or endpoints are wrong");
    require(db.sigma(excitation, 4.999) == 0.0 &&
            db.sigma(excitation, 5.0) == 0.0 &&
            near(db.sigma(excitation, 7.5), 1e-21),
            "excitation threshold or interpolation is wrong");
    require(db.sigma(ionization, 11.0) == 0.0 &&
            db.sigma(ionization, 12.0) == 0.0 &&
            near(db.sigma(ionization, 16.0), 4e-21),
            "ionization threshold or interpolation is wrong");
    require(near(db.totalSigma("N2", 0.0), 4e-20) &&
            near(db.totalSigma("N2", 10.0), 2.2e-20) &&
            near(db.totalSigma("N2", 16.0), 2.12e-20) &&
            near(db.totalSigma("O", 16.0), 1e-20),
            "totalSigma did not sum only channels of the requested species");

    const auto belowThreshold = db.activeChannels(0.0);
    const auto atTen = db.activeChannels(10.0);
    const auto atSixteen = db.activeChannels(16.0);
    require(belowThreshold.size() == 2 && atTen.size() == 3 &&
            atSixteen.size() == 4 &&
            atTen[0] == &elastic && atTen[1] == &excitation,
            "activeChannels returned wrong channels or unstable pointers");

    requireThrows<std::out_of_range>([&] {
        (void)db.sigma(elastic, 20.1);
    }, "energy above the documented range was extrapolated");
    requireThrows<std::out_of_range>([&] {
        (void)db.totalSigma("O2", 10.0);
    }, "species without a table was silently treated as zero");
    requireThrows<std::out_of_range>([&] {
        (void)db.channelById("missing");
    }, "missing channel id was accepted");
    requireThrows<std::invalid_argument>([&] {
        (void)db.activeChannels(std::numeric_limits<double>::quiet_NaN());
    }, "NaN energy was accepted");
    requireThrows<std::invalid_argument>([&] {
        (void)db.totalSigma("N2", -1.0);
    }, "negative energy was accepted");
}

void checkInvalidTables() {
    requireThrows<std::invalid_argument>([] {
        (void)TabulatedFunction({{0.0, 1e-20}});
    }, "one-point cross-section table was accepted");
    requireThrows<std::invalid_argument>([] {
        (void)TabulatedFunction({{0.0, 1e-20}, {0.0, 2e-20}});
    }, "duplicate energy was accepted");
    requireThrows<std::invalid_argument>([] {
        (void)TabulatedFunction({{0.0, 1e-20}, {10.0, -1e-20}});
    }, "negative cross section was accepted");
    requireThrows<std::invalid_argument>([] {
        (void)TabulatedFunction({{0.0, 1e-20},
            {10.0, std::numeric_limits<double>::infinity()}});
    }, "infinite cross section was accepted");
    requireThrows<std::invalid_argument>([] {
        (void)CrossSectionDatabase({});
    }, "empty cross-section database was accepted");

    auto channels = exampleChannels();
    channels[1].id = channels[0].id;
    requireThrows<std::invalid_argument>([&] {
        (void)CrossSectionDatabase(channels);
    }, "duplicate channel id was accepted");
    channels = exampleChannels();
    channels[1].source.clear();
    requireThrows<std::invalid_argument>([&] {
        (void)CrossSectionDatabase(channels);
    }, "channel without source was accepted");
    channels = exampleChannels();
    channels[1].targetState.clear();
    requireThrows<std::invalid_argument>([&] {
        (void)CrossSectionDatabase(channels);
    }, "excitation without target state was accepted");
    channels = exampleChannels();
    channels[1].threshold_eV = 6.0;
    requireThrows<std::invalid_argument>([&] {
        (void)CrossSectionDatabase(channels);
    }, "threshold not matching the first knot was accepted");
    channels = exampleChannels();
    channels[1].sigma_m2 = TabulatedFunction({{5.0, 1e-21}, {10.0, 2e-21}});
    requireThrows<std::invalid_argument>([&] {
        (void)CrossSectionDatabase(channels);
    }, "nonzero inelastic cross section at threshold was accepted");

    const CrossSectionDatabase limited({
        {"O-short", "O", CollisionType::Elastic, "", 0.0,
         TabulatedFunction({{1.0, 1e-20}, {10.0, 1e-20}}), "synthetic test"}
    });
    require(near(limited.totalSigma("O", 5.0), 1e-20),
            "elastic table with nonzero minimum energy was rejected");
    requireThrows<std::out_of_range>([&] {
        (void)limited.totalSigma("O", 0.5);
    }, "below-range elastic cross section was extrapolated");
}

void checkAtomicOxygenReference() {
    const CrossSectionDatabase db(atomicOxygenExcitationChannels());
    const auto& to1D = db.channelById("O-3P-to-1D");
    const auto& to1S = db.channelById("O-3P-to-1S");
    require(to1D.threshold_eV == 1.96 && to1S.threshold_eV == 4.18 &&
            to1D.sigma_m2.maxEnergy_eV() == 200.0 &&
            to1S.sigma_m2.maxEnergy_eV() == 200.0,
            "atomic O reference thresholds or energy limits are wrong");
    require(near(db.sigma(to1D, 5.0), 28.0e-22) &&
            near(db.sigma(to1S, 5.0), 1.60e-22) &&
            near(db.sigma(to1D, 200.0), 0.05e-22) &&
            near(db.sigma(to1S, 200.0), 0.08e-22),
            "atomic O Table 3 values or unit conversion are wrong");
    require(db.sigma(to1D, 1.95) == 0.0 &&
            db.sigma(to1S, 4.17) == 0.0 &&
            near(db.totalSigma("O", 10.0), 26.07e-22),
            "atomic O thresholds or channel sum are wrong");
    requireThrows<std::out_of_range>([&] {
        (void)db.totalSigma("O", 201.0);
    }, "atomic O data were extrapolated beyond 200 eV");
    db.requireEnergyCoverage("O", 1.0, 200.0);
    requireThrows<std::out_of_range>([&] {
        db.requireEnergyCoverage("O", 1000.0, 2000.0);
    }, "atomic O tables were claimed to cover the test electron source");
    requireThrows<std::out_of_range>([&] {
        db.requireEnergyCoverage("N2", 1.0, 10.0);
    }, "coverage accepted a species without channels");
}

void checkExplicitTails() {
    auto channels = exampleChannels();
    channels[0].tail = CrossSectionTail::InverseEnergy;
    channels[0].tailMaxEnergy_eV = 40.0;
    channels[1].tail = CrossSectionTail::Zero;
    channels[1].tailMaxEnergy_eV = 40.0;
    channels[2].tail = CrossSectionTail::Zero;
    channels[2].tailMaxEnergy_eV = 40.0;
    channels[3].tail = CrossSectionTail::Zero;
    channels[3].tailMaxEnergy_eV = 40.0;
    const CrossSectionDatabase db(channels);
    require(near(db.sigma(db.channelById("N2-elastic"), 40.0), 5e-21) &&
            db.sigma(db.channelById("N2-excitation"), 40.0) == 0.0,
            "explicit high-energy tail returned the wrong value");
    db.requireEnergyCoverage("N2", 20.0, 40.0);
    requireThrows<std::out_of_range>([&] {
        db.requireEnergyCoverage("N2", 20.0, 40.001);
    }, "tail energy limit was ignored");
    channels[0].tailMaxEnergy_eV = 20.0;
    requireThrows<std::invalid_argument>([&] {
        (void)CrossSectionDatabase(channels);
    }, "invalid model tail limit was accepted");
}

void checkLxCatImport() {
    const std::string fixture =
        "EXCITATION\n"
        "N2 -> N2(A)\n"
        " 6.0\n"
        "SPECIES: e / N2\n"
        "COLUMNS: Energy (eV) | Cross Section (m2)\n"
        "-------------------\n"
        "6.0 0\n"
        "10.0 2e-21\n"
        "-------------------\n"
        "IONIZATION\n"
        "O2 -> O2+\n"
        "12.0\n"
        "SPECIES: e / O2\n"
        "COLUMNS: Energy (eV) | Cross Section (m2)\n"
        "-------------------\n"
        "12.5 1e-22\n"
        "20.0 4e-21\n"
        "-------------------\n";
    const LxCatChannelSelection nitrogen = {
        "N2 -> N2(A)", "N2-A", "N2", CollisionType::Excitation,
        "N2(A)", "IST-Lisbon / synthetic parser fixture"};
    const LxCatChannelSelection oxygen = {
        "O2 -> O2+", "O2-ion", "O2", CollisionType::Ionization,
        "O2+", "IST-Lisbon / synthetic parser fixture"};
    std::istringstream nitrogenInput(fixture);
    std::istringstream oxygenInput(fixture);
    const CrossSectionDatabase db({
        LxCatCrossSectionLoader::load(nitrogenInput, nitrogen),
        LxCatCrossSectionLoader::load(oxygenInput, oxygen)});
    require(near(db.sigma(db.channelById("N2-A"), 8.0), 1e-21) &&
            db.sigma(db.channelById("O2-ion"), 12.0) == 0.0 &&
            near(db.sigma(db.channelById("O2-ion"), 12.5), 1e-22),
            "LXCat selection, interpolation or threshold insertion failed");
    std::istringstream missingInput(fixture);
    auto missing = nitrogen;
    missing.reaction = "N2 -> missing";
    requireThrows<std::out_of_range>([&] {
        (void)LxCatCrossSectionLoader::load(missingInput, missing);
    }, "missing LXCat process was silently accepted");
    std::istringstream badUnits(
        "EXCITATION\nN2 -> N2(A)\n6.0\n"
        "COLUMNS: Energy (eV) | Cross Section (cm2)\n"
        "-----\n6.0 0\n10.0 2e-17\n-----\n");
    requireThrows<std::invalid_argument>([&] {
        (void)LxCatCrossSectionLoader::load(badUnits, nitrogen);
    }, "LXCat wrong units were accepted");
}

} // namespace

int main() {
    try {
        checkInterpolationAndChannels();
        checkInvalidTables();
        checkAtomicOxygenReference();
        checkExplicitTails();
        checkLxCatImport();
        std::cout << "CrossSectionDatabase checks passed\n";
    } catch (const std::exception& error) {
        std::cerr << "CrossSectionDatabase check failed: " << error.what() << '\n';
        return 1;
    }
}
