#include "../src/io/IaaCrossSectionLoader.hpp"
#include "../src/atmosphere/ElectronTransportMC.hpp"
#include "../src/atmosphere/ElectronCascade.hpp"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <string>

static std::string block(const std::string& kind, const std::string& reaction,
                         const std::string& species, const std::string& parameter,
                         const std::string& units) {
    return kind + "\n" + reaction + "\n" + parameter + "\n" +
           "SPECIES: e / " + species + "\nCOLUMNS: Energy (eV) | " + units +
           "\n----------\n" +
           (kind == "ELASTIC" ? "0 2\n10 4\n" : "2 0\n10 3\n") +
           "----------\n";
}

int main() {
    std::istringstream combined(
        block("TOTAL", "O -> O", "O", "0", "Cross Section (m2)") +
        block("EFFECTIVE", "N2 -> N2", "N2", "0", "Cross Section (m2)") +
        block("EXCITATION", "O -> O(1D)", "O", "2", "Cross Section (m2)") +
        block("IONIZATION", "N2 -> N2+", "N2", "2", "Cross Section (m2)"));
    std::istringstream oxygen(block("ELASTIC", "O -> O", "O", "0",
                                    "ICS Q0 (1e-20 m2)"));
    std::istringstream oxygenMolecule(block("ELASTIC", "O2 -> O2", "O2", "0",
                                            "Cross Section (m2)"));
    std::istringstream nitrogenMolecule(block("ELASTIC", "N2 -> N2", "N2", "0",
                                              "Cross Section (m2)"));
    const auto database = IaaCrossSectionLoader::load(
        combined, oxygen, oxygenMolecule, nitrogenMolecule);
    assert(database.allChannels().size() == 5);
    assert(std::abs(database.totalSigma("O", 0) - 2e-20) < 1e-30);
    assert(database.totalSigma("O", 2) > 0);
    assert(database.channelById("IAA:O:elastic").type == CollisionType::Elastic);

    bool missingFileRejected = false;
    try {
        (void)IaaCrossSectionLoader::load("/nonexistent-aurora-iaa-data");
    } catch (const std::runtime_error&) {
        missingFileRejected = true;
    }
    assert(missingFileRejected);

    if (const char* path = std::getenv("AURORA_IAA_DATA_ROOT")) {
        const auto real = IaaCrossSectionLoader::load(path);
        assert(real.allChannels().size() > 30);
        unsigned elasticCount = 0;
        unsigned ionizationCount = 0;
        for (const auto& channel : real.allChannels()) {
            elasticCount += channel.type == CollisionType::Elastic;
            ionizationCount += channel.type == CollisionType::Ionization;
        }
        assert(elasticCount == 3);
        assert(ionizationCount == 3);
        for (const auto& species : {"O", "O2", "N2"}) {
            assert(real.channelById(std::string("IAA:") + species + ":elastic")
                       .sigma_m2.minEnergy_eV() == 0.0);
            assert(real.channelById(std::string("IAA:") + species + ":elastic")
                       .sigma_m2.maxEnergy_eV() >= 1000.0);
        }
        bool shortTableRejected = false;
        try {
            real.requireEnergyCoverage("N2", 1000.0, 2000.0);
        } catch (const std::out_of_range&) {
            shortTableRejected = true;
        }
        assert(shortTableRejected);
        const auto modeled = IaaCrossSectionLoader::load(path, 2000.0);
        for (const auto& species : {"O", "O2", "N2"}) {
            modeled.requireEnergyCoverage(species, 1000.0, 2000.0);
            assert(std::isfinite(modeled.totalSigma(species, 2000.0)));
        }
        const AtmosphereModel gas({
            {0.0, 1e18, 1e18, 1e18, 200.0},
            {2000.0, 1e18, 1e18, 1e18, 200.0}
        }, "synthetic mixed gas with local IAA tables");
        const ElectronTransportMC transport(gas, modeled);
        Random rng(123);
        const auto track = transport.run({C::earthRadius_m + 1000.0, 0.0, 0.0},
                                         {1.0, 0.0, 0.0}, 1500.0,
                                         500.0, 10.0, 1, rng);
        assert(track.events.size() == 1);
        assert(track.events[0].position_m.isFinite());
        const ElectronCascade cascade(transport, {50.0, 1.0, 100, 100, 1.0});
        Random cascadeRng(321);
        const auto cascadeResult = cascade.run(
            {C::earthRadius_m + 1000.0, 0.0, 0.0},
            {1.0, 0.0, 0.0}, 1500.0, cascadeRng);
        assert(!cascadeResult.tracks.empty());
        assert(cascadeResult.ledger.relativeResidual() < 1e-10);
        bool aboveModelRejected = false;
        try {
            modeled.requireEnergyCoverage("N2", 1000.0, 2000.1);
        } catch (const std::out_of_range&) {
            aboveModelRejected = true;
        }
        assert(aboveModelRejected);
    }
}
