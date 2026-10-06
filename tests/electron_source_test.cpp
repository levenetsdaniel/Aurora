#include "../src/ensemble/ElectronSource.hpp"
#include "../src/ensemble/ArcElectronSource.hpp"
#include "../src/physics/DipoleField.hpp"
#include "../src/physics/UniformField.hpp"

#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

template <typename Action>
void requireInvalid(Action action, const char* message) {
    try {
        action();
    } catch (const std::invalid_argument&) {
        return;
    }
    throw std::runtime_error(message);
}

void checkSampling() {
    const UniformField field(Vector3(1e-9, 2e-9, 3e-9));
    const Vector3 position(C::earthRadius_m * 2.0, 0.0, 0.0);
    constexpr std::uint64_t count = 10000;
    constexpr double flux = 3.0;
    constexpr double area_m2 = 2.0;
    const IsotropicElectronSource source(field, position, 0.1, C::electronMass_kg,
                                         1000.0, 2000.0, flux, area_m2, count);
    const ElectronSource& interface = source;
    const double expectedRate = flux * area_m2 * 1000.0 * 2.0 * std::numbers::pi;
    require(std::abs(interface.totalPhysicalRate() / expectedRate - 1.0) < 1e-14,
            "source rate does not match the analytic flux integral");
    require(source.batchSize() == count, "source lost its batch normalization size");

    Random first(12345);
    Random repeated(12345);
    const Vector3 b = field.evaluate(position) / field.evaluate(position).norm();
    double sumWeights = 0.0;
    double sumCosines = 0.0;
    double sumEnergy_eV = 0.0;
    double sumPhase_rad = 0.0;
    for (std::uint64_t index = 0; index < count; ++index) {
        const InitialCondition sample = interface.sample(first, index + 10);
        const InitialCondition again = interface.sample(repeated, index + 10);
        require(sample.id == index + 10 && sample.id == again.id &&
                sample.energy_eV == again.energy_eV &&
                sample.pitch_rad == again.pitch_rad &&
                sample.gyroPhase_rad == again.gyroPhase_rad &&
                sample.state.p.x == again.state.p.x &&
                sample.state.p.y == again.state.p.y &&
                sample.state.p.z == again.state.p.z,
                "same random seed did not reproduce the marker");
        require(sample.state.r.x == position.x && sample.state.r.y == position.y &&
                sample.state.r.z == position.z && sample.state.time_s == 0.1,
                "sampled marker has wrong starting state");
        require(sample.energy_eV >= 1000.0 && sample.energy_eV < 2000.0 &&
                sample.pitch_rad >= 0.0 && sample.pitch_rad <= std::numbers::pi / 2.0 &&
                sample.gyroPhase_rad >= 0.0 && sample.gyroPhase_rad < 2.0 * std::numbers::pi,
                "sample left the configured distribution range");

        const double cosine = -sample.state.p.dot(b) / sample.state.p.norm();
        require(std::abs(cosine - std::cos(sample.pitch_rad)) < 1e-12,
                "pitch was not measured relative to local magnetic field");
        const double energyFromMomentum_eV = Relativity::kineticEnergyFromMomentum_J(
            sample.state.p, C::electronMass_kg) / C::elementaryCharge_C;
        require(std::abs(energyFromMomentum_eV / sample.energy_eV - 1.0) < 1e-12,
                "sampled momentum does not match sampled energy");
        sumWeights += sample.weight;
        sumCosines += cosine;
        sumEnergy_eV += sample.energy_eV;
        sumPhase_rad += sample.gyroPhase_rad;
    }

    require(std::abs(sumWeights / expectedRate - 1.0) < 1e-11,
            "marker weights do not sum to physical rate");
    require(std::abs(sumCosines / count - 0.5) < 0.02,
            "isotropic hemisphere did not have uniform cos(pitch)");
    require(std::abs(sumEnergy_eV / count - 1500.0) < 15.0,
            "uniform energy distribution has the wrong mean");
    require(std::abs(sumPhase_rad / count - std::numbers::pi) < 0.06,
            "uniform gyro phase has the wrong mean");

    const IsotropicElectronSource positive(field, position, 0.0, C::electronMass_kg,
                                            1000.0, 2000.0, flux, area_m2, count, true);
    Random positiveRng(42);
    require(positive.sample(positiveRng, 1).state.p.dot(b) > 0.0,
            "positive-field hemisphere was not selected");
}

void checkInvalidParameters() {
    const UniformField field(Vector3(1e-9, 0.0, 0.0));
    const UniformField zeroField(Vector3{});
    const Vector3 position(C::earthRadius_m * 2.0, 0.0, 0.0);
    requireInvalid([&] {
        (void)IsotropicElectronSource(field, position, 0.0, C::electronMass_kg,
                                      1000.0, 2000.0, 1.0, 1.0, 0);
    }, "zero marker count was accepted");
    requireInvalid([&] {
        (void)IsotropicElectronSource(field, position, 0.0, C::electronMass_kg,
                                      2000.0, 1000.0, 1.0, 1.0, 10);
    }, "reversed energy range was accepted");
    requireInvalid([&] {
        (void)IsotropicElectronSource(field, position, 0.0, C::electronMass_kg,
                                      1000.0, 2000.0, -1.0, 1.0, 10);
    }, "negative flux was accepted");
    requireInvalid([&] {
        (void)IsotropicElectronSource(zeroField, position, 0.0, C::electronMass_kg,
                                      1000.0, 2000.0, 1.0, 1.0, 10);
    }, "zero magnetic field was accepted");
}

void checkArcSource() {
    const DipoleField field;
    constexpr std::size_t count = 10000;
    const ArcElectronSource source(field, true, 0.0, std::numbers::pi / 2.0,
                                    1e20, count);
    Random first(20260928);
    Random repeated(20260928);
    double cosineSum = 0.0;
    double weightSum = 0.0;
    for (std::size_t index = 0; index < count; ++index) {
        const auto particle = source.sample(first, index + 1);
        const auto again = source.sample(repeated, index + 1);
        const Vector3 B = field.evaluate(particle.state.r);
        const double cosine = -particle.state.p.dot(B) /
            (particle.state.p.norm() * B.norm());
        require(particle.pitch_rad == again.pitch_rad &&
                particle.state.p.x == again.state.p.x &&
                particle.state.p.y == again.state.p.y &&
                particle.state.p.z == again.state.p.z &&
                std::abs(particle.state.r.norm() - C::earthRadius_m - 450000.0) < 1e-6 &&
                particle.energy_eV >= 5000.0 && particle.energy_eV < 25000.0 &&
                particle.pitch_rad >= 0.0 &&
                particle.pitch_rad <= std::numbers::pi / 2.0 &&
                std::abs(cosine - std::cos(particle.pitch_rad)) < 1e-12,
                "arc source sampled an invalid or nonreproducible electron");
        cosineSum += cosine;
        weightSum += particle.weight;
    }
    require(std::abs(cosineSum / count - 0.5) < 0.02 &&
            std::abs(weightSum / source.totalPhysicalRate() - 1.0) < 1e-11,
            "arc source did not cover the hemisphere with normalized weights");
    requireInvalid([&] {
        (void)ArcElectronSource(field, true, 0.9, 0.8, 1e20, count);
    }, "reversed pitch range was accepted");
    requireInvalid([&] {
        (void)ArcElectronSource(field, true, 0.0, 1.6, 1e20, count);
    }, "pitch angles outside the downward hemisphere were accepted");
    const ArcElectronSource limited(field, true, 50.0 * std::numbers::pi / 180.0,
                                     80.0 * std::numbers::pi / 180.0, 1e20, count);
    Random limitedRandom(7);
    for (std::size_t index = 0; index < 100; ++index) {
        const double pitch = limited.sample(limitedRandom, index + 1).pitch_rad;
        require(pitch >= 50.0 * std::numbers::pi / 180.0 &&
                pitch <= 80.0 * std::numbers::pi / 180.0,
                "arc source ignored configured pitch limits");
    }
}

} // namespace

int main() {
    try {
        checkSampling();
        checkInvalidParameters();
        checkArcSource();
        std::cout << "ElectronSource checks passed\n";
    } catch (const std::exception& error) {
        std::cerr << "ElectronSource check failed: " << error.what() << '\n';
        return 1;
    }
}
