#ifndef AURORA_ELECTRONTRANSPORTMC_HPP
#define AURORA_ELECTRONTRANSPORTMC_HPP

#include "AtmosphereModel.hpp"
#include "CrossSectionDatabase.hpp"
#include "../math/Vector3.hpp"
#include "../numerics/Random.hpp"
#include "../physics/Constants.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

struct ElectronCollisionEvent {
    Vector3 position_m;
    double distance_m{};
    std::string channelId;
    double energyBefore_eV{};
    double energyAfter_eV{};
    CollisionType type{CollisionType::Excitation};
    Vector3 directionBefore;
    Vector3 directionAfter;
    double secondaryEnergy_eV{};
    Vector3 secondaryDirection;
};

struct ElectronTransportResult {
    std::vector<ElectronCollisionEvent> events;
    double distance_m{};
    double finalEnergy_eV{};
    bool eventLimitReached{};
    Vector3 finalPosition_m;
    Vector3 finalDirection;
    bool exitedAtmosphere{};
};

// Piecewise-straight transport with a deliberately approximate collision
// model: elastic/ionization directions are isotropic in the lab frame, elastic
// recoil is neglected, and ionization excess energy is split uniformly.
// Secondaries are returned as events, not transported (stage 18).
class ElectronTransportMC {
private:
    const AtmosphereModel& atmosphere;
    const CrossSectionDatabase& crossSections;

    static double density(const AtmosphereSample& sample,
                          const std::string& species) {
        if (species == "O") return sample.nO_m3;
        if (species == "O2") return sample.nO2_m3;
        if (species == "N2") return sample.nN2_m3;
        throw std::invalid_argument("ElectronTransportMC: unsupported species " + species);
    }

    AtmosphereSample sampleAt(const Vector3& position_m) const {
        const double altitude_m = position_m.norm() - C::earthRadius_m;
        return atmosphere.sample(altitude_m, 0.0, 0.0, 0.0);
    }

    struct Rates {
        std::vector<double> perChannel_m1;
        double total_m1{};
    };

    Rates ratesAt(const AtmosphereSample& sample, double energy_eV) const {
        for (const auto& [species, numberDensity] :
             {std::pair{"O", sample.nO_m3},
              std::pair{"O2", sample.nO2_m3},
              std::pair{"N2", sample.nN2_m3}}) {
            if (numberDensity > 0.0) {
                crossSections.requireEnergyCoverage(species, energy_eV, energy_eV);
            }
        }

        Rates rates;
        rates.perChannel_m1.reserve(crossSections.allChannels().size());
        for (const CollisionChannel& channel : crossSections.allChannels()) {
            const double numberDensity = density(sample, channel.species);
            const double rate_m1 = numberDensity == 0.0 ? 0.0 :
                numberDensity * crossSections.sigma(channel, energy_eV);
            if (!std::isfinite(rate_m1) || rate_m1 < 0.0) {
                throw std::overflow_error("ElectronTransportMC: collision rate overflow");
            }
            rates.perChannel_m1.push_back(rate_m1);
            rates.total_m1 += rate_m1;
        }
        if (!std::isfinite(rates.total_m1)) {
            throw std::overflow_error("ElectronTransportMC: total collision rate overflow");
        }
        return rates;
    }

    static double sampleOpticalDepth(Random& rng) {
        double uniform;
        do {
            uniform = rng.uniform01();
        } while (uniform == 0.0);
        return -std::log(uniform);
    }

    static Vector3 isotropicDirection(Random& rng) {
        const double cosine = 2.0 * rng.uniform01() - 1.0;
        const double azimuth = 2.0 * std::numbers::pi * rng.uniform01();
        const double sine = std::sqrt(std::max(0.0, 1.0 - cosine * cosine));
        return {sine * std::cos(azimuth), sine * std::sin(azimuth), cosine};
    }

    double boundaryDistance(const Vector3& position_m, const Vector3& direction) const {
        const double radius_m = position_m.norm();
        const double inner_m = C::earthRadius_m + atmosphere.minAltitude_m();
        const double outer_m = C::earthRadius_m + atmosphere.maxAltitude_m();
        const double radialDirection = position_m.dot(direction);
        constexpr double surfaceTolerance_m = 1e-6;
        if ((std::abs(radius_m - inner_m) <= surfaceTolerance_m && radialDirection < 0.0) ||
            (std::abs(radius_m - outer_m) <= surfaceTolerance_m && radialDirection > 0.0)) {
            return 0.0;
        }
        double nearest_m = std::numeric_limits<double>::infinity();
        for (const double surface_m : {inner_m, outer_m}) {
            const double discriminant = radialDirection * radialDirection -
                (position_m.norm2() - surface_m * surface_m);
            if (discriminant < 0.0) continue;
            const double root = std::sqrt(discriminant);
            for (const double distance_m : {-radialDirection - root,
                                            -radialDirection + root}) {
                if (distance_m > surfaceTolerance_m) {
                    nearest_m = std::min(nearest_m, distance_m);
                }
            }
        }
        return nearest_m;
    }

public:
    ElectronTransportMC(const AtmosphereModel& atmosphere,
                        const CrossSectionDatabase& crossSections)
        : atmosphere(atmosphere), crossSections(crossSections) {
        for (const CollisionChannel& channel : crossSections.allChannels()) {
            if (channel.species != "O" && channel.species != "O2" &&
                channel.species != "N2") {
                throw std::invalid_argument("ElectronTransportMC: unsupported species");
            }
        }
    }

    ElectronTransportResult run(const Vector3& start_m, const Vector3& direction,
                                double initialEnergy_eV, double pathLength_m,
                                double integrationStep_m, std::size_t maxEvents,
                                Random& rng) const {
        if (!start_m.isFinite() || !direction.isFinite() ||
            !std::isfinite(direction.norm()) || direction.norm() <= 0.0 ||
            !std::isfinite(initialEnergy_eV) || initialEnergy_eV < 0.0 ||
            !std::isfinite(pathLength_m) || pathLength_m < 0.0 ||
            !std::isfinite(integrationStep_m) || integrationStep_m <= 0.0 ||
            maxEvents == 0) {
            throw std::invalid_argument("ElectronTransportMC: invalid track parameters");
        }
        const Vector3 unitDirection = direction / direction.norm();
        (void)sampleAt(start_m);

        ElectronTransportResult result{{}, 0.0, initialEnergy_eV, false,
                                       start_m, unitDirection, false};
        double remainingOpticalDepth = sampleOpticalDepth(rng);
        while (result.distance_m < pathLength_m) {
            const double boundary_m = boundaryDistance(result.finalPosition_m,
                                                       result.finalDirection);
            if (boundary_m == 0.0) {
                result.exitedAtmosphere = true;
                return result;
            }
            const double step_m = std::min({integrationStep_m,
                                            pathLength_m - result.distance_m,
                                            boundary_m});
            const Vector3 midpoint = result.finalPosition_m +
                result.finalDirection * (step_m / 2.0);
            const Rates rates = ratesAt(sampleAt(midpoint), result.finalEnergy_eV);
            const double stepOpticalDepth = rates.total_m1 * step_m;
            if (!std::isfinite(stepOpticalDepth)) {
                throw std::overflow_error("ElectronTransportMC: optical depth overflow");
            }
            if (rates.total_m1 == 0.0 || stepOpticalDepth < remainingOpticalDepth) {
                remainingOpticalDepth -= stepOpticalDepth;
                result.finalPosition_m = result.finalPosition_m + result.finalDirection * step_m;
                result.distance_m += step_m;
                if (step_m == boundary_m) {
                    result.exitedAtmosphere = true;
                    return result;
                }
                continue;
            }

            const double collisionDistance_m = remainingOpticalDepth / rates.total_m1;
            result.finalPosition_m = result.finalPosition_m +
                result.finalDirection * collisionDistance_m;
            result.distance_m += collisionDistance_m;
            const double draw_m1 = rng.uniform01() * rates.total_m1;
            double cumulative_m1 = 0.0;
            std::size_t selected = rates.perChannel_m1.size();
            for (std::size_t index = 0; index < rates.perChannel_m1.size(); ++index) {
                cumulative_m1 += rates.perChannel_m1[index];
                if (draw_m1 < cumulative_m1) {
                    selected = index;
                    break;
                }
            }
            if (selected == rates.perChannel_m1.size()) {
                throw std::runtime_error("ElectronTransportMC: failed to select collision channel");
            }
            const CollisionChannel& channel = crossSections.allChannels()[selected];
            const double energyBefore_eV = result.finalEnergy_eV;
            const Vector3 directionBefore = result.finalDirection;
            double secondaryEnergy_eV = 0.0;
            Vector3 secondaryDirection;
            switch (channel.type) {
                case CollisionType::Elastic:
                    result.finalDirection = isotropicDirection(rng);
                    break;
                case CollisionType::Excitation:
                    result.finalEnergy_eV = std::max(0.0,
                        energyBefore_eV - channel.threshold_eV);
                    break;
                case CollisionType::Ionization: {
                    const double excess_eV = std::max(0.0,
                        energyBefore_eV - channel.threshold_eV);
                    result.finalEnergy_eV = excess_eV * rng.uniform01();
                    secondaryEnergy_eV = excess_eV - result.finalEnergy_eV;
                    result.finalDirection = isotropicDirection(rng);
                    secondaryDirection = isotropicDirection(rng);
                    break;
                }
            }
            result.events.push_back({result.finalPosition_m, result.distance_m, channel.id,
                                     energyBefore_eV, result.finalEnergy_eV, channel.type,
                                     directionBefore, result.finalDirection,
                                     secondaryEnergy_eV, secondaryDirection});
            if (result.events.size() == maxEvents) {
                result.eventLimitReached = result.distance_m < pathLength_m &&
                    boundaryDistance(result.finalPosition_m, result.finalDirection) > 0.0;
                return result;
            }
            remainingOpticalDepth = sampleOpticalDepth(rng);
        }
        return result;
    }
};

#endif // AURORA_ELECTRONTRANSPORTMC_HPP
