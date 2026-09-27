#include "../src/atmosphere/AtmosphereModel.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

bool near(double actual, double expected, double relativeTolerance = 1e-12) {
    return std::abs(actual - expected) <= relativeTolerance * std::abs(expected);
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

std::vector<AtmosphereSample> table() {
    return {
        {100000.0, 1e18, 4e18, 9e18, 200.0},
        {120000.0, 1e16, 4e16, 9e16, 400.0},
        {140000.0, 1e14, 4e14, 9e14, 600.0}
    };
}

void checkInterpolation() {
    const AtmosphereModel atmosphere(table(), "synthetic test profile, SI units");
    require(atmosphere.minAltitude_m() == 100000.0 &&
            atmosphere.maxAltitude_m() == 140000.0 &&
            atmosphere.source() == "synthetic test profile, SI units",
            "model lost altitude range or source metadata");

    const auto lower = atmosphere.sample(100000.0, 0.0, 0.0, 0.0);
    const auto node = atmosphere.sample(120000.0, 0.1, 1.0, 10.0);
    const auto upper = atmosphere.sample(140000.0, -0.1, -1.0, 20.0);
    require(lower.nO_m3 == 1e18 && lower.nO2_m3 == 4e18 &&
            lower.nN2_m3 == 9e18 && lower.temperature_K == 200.0 &&
            node.nO_m3 == 1e16 && node.temperature_K == 400.0 &&
            upper.nO_m3 == 1e14 && upper.temperature_K == 600.0,
            "table nodes were not preserved exactly");

    const auto midpoint = atmosphere.sample(110000.0, 0.0, 0.0, 0.0);
    require(midpoint.altitude_m == 110000.0 &&
            near(midpoint.nO_m3, 1e17) &&
            near(midpoint.nO2_m3, 4e17) &&
            near(midpoint.nN2_m3, 9e17) &&
            midpoint.temperature_K == 300.0,
            "log-density or linear-temperature interpolation is wrong");
    require(!near(midpoint.nO_m3, (1e18 + 1e16) / 2.0, 1e-6),
            "density was interpolated linearly instead of logarithmically");

    const auto otherCoordinates = atmosphere.sample(110000.0,
        std::numbers::pi / 3.0, -2.0, 1000.0);
    require(otherCoordinates.nO_m3 == midpoint.nO_m3 &&
            otherCoordinates.nO2_m3 == midpoint.nO2_m3 &&
            otherCoordinates.nN2_m3 == midpoint.nN2_m3 &&
            otherCoordinates.temperature_K == midpoint.temperature_K,
            "1D atmosphere incorrectly depends on latitude, longitude, or time");

    const auto justBelow = atmosphere.sample(120000.0 - 1e-6, 0.0, 0.0, 0.0);
    const auto justAbove = atmosphere.sample(120000.0 + 1e-6, 0.0, 0.0, 0.0);
    require(justBelow.nO_m3 > 0.0 && justAbove.nO_m3 > 0.0 &&
            near(justBelow.nO_m3, node.nO_m3, 1e-8) &&
            near(justAbove.nO_m3, node.nO_m3, 1e-8),
            "density is negative or discontinuous near a table node");

    const AtmosphereModel extreme({
        {0.0, 1e-300, 1e-300, 1e-300, 100.0},
        {1000.0, 1e300, 1e300, 1e300, 200.0}
    }, "synthetic numerical-range test");
    require(near(extreme.sample(500.0, 0.0, 0.0, 0.0).nO_m3, 1.0),
            "log interpolation is unstable across a large density range");

    const AtmosphereModel zeroSpecies({
        {0.0, 0.0, 1e10, 0.0, 100.0},
        {1000.0, 2e10, 1e8, 0.0, 200.0}
    }, "synthetic zero-density test");
    const auto zeroMidpoint = zeroSpecies.sample(500.0, 0.0, 0.0, 0.0);
    require(near(zeroMidpoint.nO_m3, 1e10) &&
            near(zeroMidpoint.nO2_m3, 1e9) &&
            zeroMidpoint.nN2_m3 == 0.0,
            "zero-density intervals were not interpolated safely");
}

void checkInvalidInputs() {
    const auto valid = table();
    requireThrows<std::invalid_argument>([] {
        (void)AtmosphereModel({}, "empty");
    }, "empty atmosphere table was accepted");
    requireThrows<std::invalid_argument>([&] {
        (void)AtmosphereModel({valid.front()}, "one point");
    }, "one-point atmosphere table was accepted");
    requireThrows<std::invalid_argument>([&] {
        (void)AtmosphereModel(valid, "");
    }, "missing source metadata was accepted");

    auto invalid = valid;
    invalid[1].altitude_m = invalid[0].altitude_m;
    requireThrows<std::invalid_argument>([&] {
        (void)AtmosphereModel(invalid, "duplicate altitude");
    }, "duplicate altitude was accepted");
    invalid = valid;
    invalid[0].altitude_m = -1.0;
    requireThrows<std::invalid_argument>([&] {
        (void)AtmosphereModel(invalid, "negative altitude");
    }, "negative altitude was accepted");
    invalid = valid;
    invalid[1].nO2_m3 = -1.0;
    requireThrows<std::invalid_argument>([&] {
        (void)AtmosphereModel(invalid, "negative density");
    }, "negative density was accepted");
    invalid = valid;
    invalid[1].temperature_K = std::numeric_limits<double>::infinity();
    requireThrows<std::invalid_argument>([&] {
        (void)AtmosphereModel(invalid, "infinite temperature");
    }, "infinite temperature was accepted");

    const AtmosphereModel atmosphere(valid, "synthetic test");
    requireThrows<std::out_of_range>([&] {
        (void)atmosphere.sample(99999.0, 0.0, 0.0, 0.0);
    }, "altitude below table was extrapolated");
    requireThrows<std::out_of_range>([&] {
        (void)atmosphere.sample(140001.0, 0.0, 0.0, 0.0);
    }, "altitude above table was extrapolated");
    requireThrows<std::invalid_argument>([&] {
        (void)atmosphere.sample(std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0, 0.0);
    }, "NaN altitude was accepted");
    requireThrows<std::invalid_argument>([&] {
        (void)atmosphere.sample(110000.0, std::numbers::pi, 0.0, 0.0);
    }, "invalid latitude was accepted");
    requireThrows<std::invalid_argument>([&] {
        (void)atmosphere.sample(110000.0, 0.0, 0.0, -1.0);
    }, "negative time was accepted");
}

} // namespace

int main() {
    try {
        checkInterpolation();
        checkInvalidInputs();
        std::cout << "AtmosphereModel checks passed\n";
    } catch (const std::exception& error) {
        std::cerr << "AtmosphereModel check failed: " << error.what() << '\n';
        return 1;
    }
}
