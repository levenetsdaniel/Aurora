#include "../src/io/AtmosphereCsvLoader.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

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

void checkReferenceProfile() {
    const std::filesystem::path path =
        std::filesystem::path(AURORA_SOURCE_DIR) / "data/atmosphere" /
        "nrlmsis20_2020-03-20_70N_0E_demo.csv";
    const AtmosphereModel model = AtmosphereCsvLoader::load(path);
    require(model.minAltitude_m() == 70000.0 && model.maxAltitude_m() == 300000.0,
            "reference profile has wrong altitude range");
    require(model.source().find("NRLMSIS 2.0") != std::string::npos &&
            model.source().find("f107_daily_previous_day_sfu=70") != std::string::npos &&
            model.source().find("geodetic_latitude_deg=70") != std::string::npos,
            "reference profile lost scenario provenance");

    const AtmosphereSample node = model.sample(150000.0, 0.0, 0.0, 0.0);
    require(std::abs(node.nO_m3 / 1.0896557160949110e16 - 1.0) < 1e-13 &&
            std::abs(node.nO2_m3 / 3.0736153792486445e15 - 1.0) < 1e-13 &&
            std::abs(node.nN2_m3 / 2.1740890998269836e16 - 1.0) < 1e-13 &&
            std::abs(node.temperature_K - 627.3459771724722) < 1e-10,
            "reference profile SI columns were read incorrectly");

    const AtmosphereSample upper = model.sample(155000.0, 0.0, 0.0, 0.0);
    const AtmosphereSample midpoint = model.sample(152500.0, 0.0, 0.0, 0.0);
    require(std::abs(midpoint.nO_m3 /
            std::sqrt(node.nO_m3 * upper.nO_m3) - 1.0) < 1e-12 &&
            std::abs(midpoint.temperature_K -
            (node.temperature_K + upper.temperature_K) * 0.5) < 1e-12,
            "loaded profile did not interpolate between table nodes");
}

void checkInvalidCsv() {
    constexpr const char* header =
        "altitude_m,nO_m3,nO2_m3,nN2_m3,temperature_K\n";
    const auto parse = [](const std::string& contents) {
        std::istringstream input(contents);
        return AtmosphereCsvLoader::parse(input);
    };
    const std::string valid = std::string("# source=synthetic test\n") + header +
        "100000,1e16,2e16,3e16,200\n" +
        "110000,1e15,2e15,3e15,300\n";
    require(parse(valid).sample(105000.0, 0.0, 0.0, 0.0).nO_m3 > 0.0,
            "valid in-memory CSV was rejected");

    requireThrows<std::runtime_error>([&] {
        (void)parse(std::string(header) +
                    "100000,1e16,2e16,3e16,200\n110000,1e15,2e15,3e15,300\n");
    }, "missing source metadata was accepted");
    requireThrows<std::runtime_error>([&] {
        (void)parse("# source=synthetic test\naltitude_km,nO_cm3,nO2_cm3,nN2_cm3,temperature_K\n"
                    "100,1e10,2e10,3e10,200\n110,1e9,2e9,3e9,300\n");
    }, "wrong-unit CSV header was accepted");
    requireThrows<std::runtime_error>([&] {
        (void)parse(std::string("# source=synthetic test\n") + header +
                    "100000,1e16,2e16,3e16,200\n110000,1e15,2e15,300\n");
    }, "missing CSV column was accepted");
    requireThrows<std::runtime_error>([&] {
        (void)parse(std::string("# source=synthetic test\n") + header +
                    "100000,1e16,2e16,3e16,200,extra\n110000,1e15,2e15,3e15,300\n");
    }, "extra CSV column was accepted");
    requireThrows<std::runtime_error>([&] {
        (void)parse(std::string("# source=synthetic test\n") + header +
                    "100000,1e16,2e16,3e16,200\n110000,invalid,2e15,3e15,300\n");
    }, "non-numeric CSV value was accepted");
    requireThrows<std::invalid_argument>([&] {
        (void)parse(std::string("# source=synthetic test\n") + header +
                    "110000,1e16,2e16,3e16,200\n100000,1e15,2e15,3e15,300\n");
    }, "unordered altitudes were accepted");
    requireThrows<std::runtime_error>([] {
        (void)AtmosphereCsvLoader::load("/private/tmp/aurora-nonexistent-profile.csv");
    }, "missing atmosphere file was accepted");
}

} // namespace

int main() {
    try {
        checkReferenceProfile();
        checkInvalidCsv();
        std::cout << "AtmosphereCsvLoader checks passed\n";
    } catch (const std::exception& error) {
        std::cerr << "AtmosphereCsvLoader check failed: " << error.what() << '\n';
        return 1;
    }
}
