#include "../src/visualization/AuroraPipeline.hpp"
#include "../src/visualization/ProgressReporter.hpp"

#include <atomic>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

CrossSectionDatabase crossSections() {
    return CrossSectionDatabase({
        {"red", "O", CollisionType::Excitation, "O(1D)", 2.0,
         TabulatedFunction({{2.0, 0.0}, {3.0, 1e-18}, {20.0, 1e-18}}), "test"},
        {"green", "O", CollisionType::Excitation, "O(1S)", 4.0,
         TabulatedFunction({{4.0, 0.0}, {5.0, 1e-18}, {20.0, 1e-18}}), "test"},
        {"other", "O", CollisionType::Excitation, "unmodelled", 2.0,
         TabulatedFunction({{2.0, 0.0}, {3.0, 5e-19}, {20.0, 5e-19}}), "test"}
    });
}

KineticsCatalog kineticsCatalog() {
    return KineticsCatalog(
        {{"O(1D)", "O", 2.0, "test"}, {"O(1S)", "O", 4.0, "test"}},
        {{"O(1D)", "O", 0, 300, 0, 0, "test"},
         {"O(1D)", "O2", 0, 300, 0, 0, "test"},
         {"O(1D)", "N2", 0, 300, 0, 0, "test"},
         {"O(1S)", "O", 0, 300, 0, 0, "test"},
         {"O(1S)", "O2", 0, 300, 0, 0, "test"},
         {"O(1S)", "N2", 0, 300, 0, 0, "test"}},
        {{{"O", "O(1D)"}, "O(1D)"}, {{"O", "O(1S)"}, "O(1S)"}});
}

void checkParallelReproducibility() {
    const AtmosphereModel atmosphere({
        {70000, 1e18, 0, 0, 300}, {300000, 1e18, 0, 0, 300}}, "test");
    const auto cross = crossSections();
    const auto states = kineticsCatalog();
    const SpectralLineCatalog lines({
        {"O_5577", 557.7, "O(1S)", 3.0, "test"},
        {"O_6300", 630.0, "O(1D)", 1.0, "test"}}, states);
    const double radius_m = C::earthRadius_m + 150000;
    std::vector<PrecipitatingElectron> electrons;
    for (std::uint64_t id = 1; id <= 8; ++id) {
        const double momentum = Relativity::momentumMagnitudeFromKineticEnergy_J(
            10.0 * C::elementaryCharge_C, C::electronMass_kg);
        electrons.push_back({id, {{radius_m, 0, 0}, {-momentum, 0, 0}, 0},
                             10.0, 1.0, 5.0});
    }
    AuroraPipelineConfig config{
        {{radius_m - 300, -500, -500}, {100, 1000, 1000}, 4, 1, 1},
        {100, 1, 100, 10, 1}, 1, 3, 12345};
    const auto serial = AuroraPipeline(atmosphere, cross, states, lines, config).run(electrons);
    config.threads = 4;
    config.batchSize = 8;
    std::atomic<std::size_t> completed{0};
    const auto parallel = AuroraPipeline(atmosphere, cross, states, lines, config).run(
        electrons, [&] { completed.fetch_add(1, std::memory_order_relaxed); });
    require(serial.primaryCount == electrons.size() &&
            completed.load(std::memory_order_relaxed) == electrons.size() &&
            parallel.primaryCount == serial.primaryCount &&
            parallel.cascadeTracks == serial.cascadeTracks &&
            parallel.excitationEvents == serial.excitationEvents &&
            parallel.excitationEvents > 0 &&
            parallel.initialEnergy_eV_s1 == serial.initialEnergy_eV_s1 &&
            std::abs(parallel.energyResidual_eV_s1) < 1e-8 &&
            serial.untrackedExcitationRate_s1 > 0.0 &&
            serial.untrackedExcitationRate_s1 ==
                parallel.untrackedExcitationRate_s1 &&
            serial.outsideGridRate_s1 == 0.0,
            "parallel cascade diagnostics differ from serial run");
    double total = 0.0;
    for (std::size_t x = 0; x < 4; ++x) {
        const double red = serial.emission.emissivity("O_6300", x, 0, 0);
        const double green = serial.emission.emissivity("O_5577", x, 0, 0);
        require(red == parallel.emission.emissivity("O_6300", x, 0, 0) &&
                green == parallel.emission.emissivity("O_5577", x, 0, 0),
                "emission changed with thread count or batch size");
        total += red + green;
    }
    require(total > 0.0, "cascade produced no visible emission");
    std::ostringstream csv;
    parallel.emission.writeSparseCsv(csv);
    require(csv.str().find("j_O_5577_ph_m3_s1,j_O_6300_ph_m3_s1") !=
            std::string::npos,
            "emission CSV lacks separate physical spectral fields");
}

void checkProgressReporter() {
    std::ostringstream output;
    ProgressReporter progress("Atmospheric cascades", 20, output);
    parallelFor(20, 4, [&](std::size_t) { progress.advance(); });
    const std::string log = output.str();
    require(log.find("Atmospheric cascades: 0% (0/20)\n") == 0 &&
            log.find("Atmospheric cascades: 100% (20/20)\n") != std::string::npos,
            "progress reporter did not log start and completion");
}

} // namespace

int main() {
    try {
        checkParallelReproducibility();
        checkProgressReporter();
        std::cout << "AuroraPipeline checks passed\n";
    } catch (const std::exception& error) {
        std::cerr << "AuroraPipeline check failed: " << error.what() << '\n';
        return 1;
    }
}
