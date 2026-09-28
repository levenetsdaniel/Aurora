#include "../src/atmosphere/ElectronTransportMC.hpp"

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

AtmosphereModel atmosphere(double nO_m3) {
    return AtmosphereModel({
        {0.0, nO_m3, 0.0, 0.0, 200.0},
        {2000.0, nO_m3, 0.0, 0.0, 200.0}
    }, "synthetic homogeneous oxygen atmosphere");
}

CrossSectionDatabase channels() {
    return CrossSectionDatabase({
        {"O-A", "O", CollisionType::Excitation, "A", 1.0,
         TabulatedFunction({{1.0, 0.0}, {2.0, 1e-20}, {20.0, 1e-20}}),
         "synthetic test"},
        {"O-B", "O", CollisionType::Excitation, "B", 2.0,
         TabulatedFunction({{2.0, 0.0}, {3.0, 3e-20}, {20.0, 3e-20}}),
         "synthetic test"}
    });
}

void checkHomogeneousStatistics() {
    const auto gas = atmosphere(1e18);
    const auto crossSections = channels();
    const ElectronTransportMC transport(gas, crossSections);
    const Vector3 start(C::earthRadius_m + 1500.0, 0.0, 0.0);
    const Vector3 inward(-1.0, 0.0, 0.0);
    Random rng(123456);
    constexpr int count = 20000;
    double sumDistance_m = 0.0;
    int selectedA = 0;
    for (int index = 0; index < count; ++index) {
        const auto result = transport.run(start, inward, 10.0, 1000.0,
                                          25.0, 1, rng);
        require(result.events.size() == 1 && result.eventLimitReached &&
                result.events.front().distance_m >= 0.0 &&
                result.events.front().distance_m <= 1000.0 &&
                result.events.front().position_m.isFinite(),
                "homogeneous track did not produce a valid collision");
        sumDistance_m += result.events.front().distance_m;
        selectedA += result.events.front().channelId == "O-A";
        const double loss_eV = result.events.front().channelId == "O-A" ? 1.0 : 2.0;
        require(result.finalEnergy_eV == 10.0 - loss_eV,
                "excitation did not subtract its threshold energy");
    }
    // n * (sigma_A + sigma_B) = 0.04 m^-1, so mean free path = 25 m.
    require(std::abs(sumDistance_m / count - 25.0) < 0.6,
            "sampled mean free path differs from 1 / total rate");
    require(std::abs(static_cast<double>(selectedA) / count - 0.25) < 0.015,
            "channel frequencies differ from relative cross sections");

    int noCollision = 0;
    for (int index = 0; index < count; ++index) {
        noCollision += transport.run(start, inward, 10.0, 10.0,
                                     2.0, 1, rng).events.empty();
    }
    require(std::abs(static_cast<double>(noCollision) / count - std::exp(-0.4)) < 0.015,
            "no-collision probability differs from exp(-n sigma L)");
}

void checkZeroAndDepletedEnergy() {
    const auto emptyGas = atmosphere(0.0);
    const auto crossSections = channels();
    const ElectronTransportMC emptyTransport(emptyGas, crossSections);
    Random rng(7);
    const Vector3 start(C::earthRadius_m + 1500.0, 0.0, 0.0);
    const auto emptyResult = emptyTransport.run(start, {-1.0, 0.0, 0.0},
                                                 10.0, 1000.0, 100.0, 10, rng);
    require(emptyResult.events.empty() && emptyResult.distance_m == 1000.0 &&
            emptyResult.finalEnergy_eV == 10.0 && !emptyResult.eventLimitReached,
            "zero-density track did not escape unchanged");

    const CrossSectionDatabase oneChannel({
        {"O-loss", "O", CollisionType::Excitation, "loss", 5.0,
         TabulatedFunction({{5.0, 0.0}, {6.0, 1e-20}, {20.0, 1e-20}}),
         "synthetic test"}
    });
    const auto gas = atmosphere(1e20);
    const ElectronTransportMC transport(gas, oneChannel);
    const auto result = transport.run(start, {-1.0, 0.0, 0.0},
                                      10.0, 1000.0, 20.0, 10, rng);
    require(result.events.size() == 1 && result.finalEnergy_eV == 5.0 &&
            result.distance_m == 1000.0 && !result.eventLimitReached,
            "track did not stop colliding after reaching excitation threshold");
}

void checkElasticAndIonization() {
    const auto gas = atmosphere(1e18);
    const Vector3 start(C::earthRadius_m + 1000.0, 0.0, 0.0);
    const CrossSectionDatabase elastic({
        {"O-elastic", "O", CollisionType::Elastic, "", 0.0,
         TabulatedFunction({{0.0, 1e-20}, {20.0, 1e-20}}), "synthetic test"}
    });
    const ElectronTransportMC elasticTransport(gas, elastic);
    Random elasticRng(314);
    const auto scattered = elasticTransport.run(start, {1.0, 0.0, 0.0},
                                                 10.0, 500.0, 10.0, 1, elasticRng);
    require(scattered.events.size() == 1 &&
            scattered.events[0].type == CollisionType::Elastic &&
            scattered.finalEnergy_eV == 10.0 &&
            std::abs(scattered.finalDirection.norm() - 1.0) < 1e-12 &&
            (scattered.finalDirection - Vector3{1.0, 0.0, 0.0}).norm() > 1e-6 &&
            scattered.events[0].secondaryEnergy_eV == 0.0,
            "elastic collision did not change direction without energy loss");

    const CrossSectionDatabase ionization({
        {"O-ionization", "O", CollisionType::Ionization, "", 5.0,
         TabulatedFunction({{5.0, 0.0}, {6.0, 1e-20}, {20.0, 1e-20}}),
         "synthetic test"}
    });
    const ElectronTransportMC ionizationTransport(gas, ionization);
    Random ionizationRng(2718);
    const auto ionized = ionizationTransport.run(start, {1.0, 0.0, 0.0},
                                                 10.0, 500.0, 10.0, 1,
                                                 ionizationRng);
    require(ionized.events.size() == 1 &&
            ionized.events[0].type == CollisionType::Ionization &&
            ionized.events[0].secondaryEnergy_eV > 0.0 &&
            std::abs(ionized.finalEnergy_eV +
                     ionized.events[0].secondaryEnergy_eV + 5.0 - 10.0) < 1e-12 &&
            std::abs(ionized.events[0].secondaryDirection.norm() - 1.0) < 1e-12,
            "ionization did not split excess energy into primary and secondary");
}

void checkScatteringAndBoundary() {
    const auto gas = atmosphere(1e18);
    const CrossSectionDatabase elastic({
        {"O-elastic", "O", CollisionType::Elastic, "", 0.0,
         TabulatedFunction({{0.0, 1e-20}, {20.0, 1e-20}}), "synthetic test"}
    });
    const ElectronTransportMC transport(gas, elastic);
    const Vector3 start(C::earthRadius_m + 1000.0, 0.0, 0.0);
    Random rng(42);
    const auto result = transport.run(start, {1.0, 0.0, 0.0},
                                      10.0, 200.0, 5.0, 10, rng);
    require(result.events.size() >= 2 && result.finalPosition_m.isFinite() &&
            std::abs(result.finalDirection.norm() - 1.0) < 1e-12,
            "track did not continue along scattered directions");
    for (std::size_t index = 1; index < result.events.size(); ++index) {
        const auto& previous = result.events[index - 1];
        const auto& next = result.events[index];
        const Vector3 displacement = next.position_m - previous.position_m;
        require(std::abs(displacement.norm() -
                         (next.distance_m - previous.distance_m)) < 1e-6 &&
                displacement.dot(previous.directionAfter) > 0.0,
                "post-collision track did not follow the new direction");
    }
    Random exitRng(2);
    const auto exiting = transport.run({C::earthRadius_m + 1990.0, 0.0, 0.0},
                                       {1.0, 0.0, 0.0}, 10.0, 100.0,
                                       20.0, 10, exitRng);
    require(exiting.exitedAtmosphere && exiting.distance_m <= 10.001 &&
            exiting.finalPosition_m.isFinite(),
            "track crossing atmosphere boundary was not stopped");
}

void checkMonteCarloConvergence() {
    const auto gas = atmosphere(1e18);
    const auto crossSections = channels();
    const ElectronTransportMC transport(gas, crossSections);
    const Vector3 start(C::earthRadius_m + 1500.0, 0.0, 0.0);
    const double expected = std::exp(-0.4);
    auto rootMeanSquareError = [&](int tracksPerBatch) {
        double squaredError = 0.0;
        for (int batch = 0; batch < 8; ++batch) {
            Random rng(9000 + batch);
            int noCollision = 0;
            for (int track = 0; track < tracksPerBatch; ++track) {
                noCollision += transport.run(start, {-1.0, 0.0, 0.0},
                                             10.0, 10.0, 2.0, 1, rng).events.empty();
            }
            const double error = static_cast<double>(noCollision) / tracksPerBatch - expected;
            squaredError += error * error;
        }
        return std::sqrt(squaredError / 8.0);
    };
    require(rootMeanSquareError(5000) < rootMeanSquareError(500),
            "Monte Carlo probability estimate did not converge with more tracks");
}

void checkRejections() {
    const auto gas = atmosphere(1e18);
    const auto crossSections = channels();
    const ElectronTransportMC transport(gas, crossSections);
    const Vector3 start(C::earthRadius_m + 1500.0, 0.0, 0.0);
    Random rng(42);
    requireThrows<std::out_of_range>([&] {
        (void)transport.run(start, {-1.0, 0.0, 0.0},
                            1000.0, 100.0, 10.0, 1, rng);
    }, "energy beyond table was silently extrapolated");
    requireThrows<std::invalid_argument>([&] {
        (void)transport.run(start, {}, 10.0, 100.0, 10.0, 1, rng);
    }, "zero direction was accepted");
    requireThrows<std::invalid_argument>([&] {
        (void)transport.run(start, {-1.0, 0.0, 0.0},
                            10.0, 100.0, 0.0, 1, rng);
    }, "zero integration step was accepted");
    const auto molecularGas = AtmosphereModel({
        {0.0, 1e18, 1e18, 0.0, 200.0},
        {2000.0, 1e18, 1e18, 0.0, 200.0}
    }, "synthetic mixed atmosphere");
    const ElectronTransportMC incomplete(molecularGas, crossSections);
    requireThrows<std::out_of_range>([&] {
        (void)incomplete.run(start, {-1.0, 0.0, 0.0},
                             10.0, 100.0, 10.0, 1, rng);
    }, "gas with no O2 channels was silently ignored");
}

} // namespace

int main() {
    try {
        checkHomogeneousStatistics();
        checkZeroAndDepletedEnergy();
        checkElasticAndIonization();
        checkScatteringAndBoundary();
        checkMonteCarloConvergence();
        checkRejections();
        std::cout << "ElectronTransportMC checks passed\n";
    } catch (const std::exception& error) {
        std::cerr << "ElectronTransportMC check failed: " << error.what() << '\n';
        return 1;
    }
}
