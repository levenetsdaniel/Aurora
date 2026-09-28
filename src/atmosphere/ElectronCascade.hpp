#ifndef AURORA_ELECTRONCASCADE_HPP
#define AURORA_ELECTRONCASCADE_HPP

#include "ElectronTransportMC.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <deque>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

struct EnergyLedger {
    double initial_eV{};
    double excitation_eV{};
    double ionization_eV{};
    double thermal_eV{};
    double escaped_eV{};
    double remaining_eV{};
    double secondaryTransfer_eV{};

    double residual_eV() const {
        return initial_eV - excitation_eV - ionization_eV - thermal_eV -
               escaped_eV - remaining_eV - secondaryTransfer_eV;
    }

    double relativeResidual() const {
        return std::abs(residual_eV()) / std::max(1e-12, initial_eV);
    }
};

struct ElectronCascadeConfig {
    double pathLength_m{};
    double integrationStep_m{};
    std::size_t maxEventsPerTrack{};
    std::size_t maxTracks{};
    double thermalizationEnergy_eV{};
};

struct ElectronTrackRecord {
    static constexpr std::size_t noParent = std::numeric_limits<std::size_t>::max();

    std::size_t parentIndex{noParent};
    std::size_t generation{};
    double initialEnergy_eV{};
    ElectronTransportResult transport;
    EnergyLedger ledger;
    bool thermalizedWithoutTransport{};
};

struct ElectronCascadeResult {
    std::vector<ElectronTrackRecord> tracks;
    EnergyLedger ledger;
    double secondaryGenerated_eV{};
    std::size_t unprocessedTracks{};
    bool limitReached{};
};

class ElectronCascade {
private:
    struct PendingTrack {
        Vector3 position_m;
        Vector3 direction;
        double energy_eV{};
        std::size_t parentIndex{ElectronTrackRecord::noParent};
        std::size_t generation{};
    };

    const ElectronTransportMC& transport;
    ElectronCascadeConfig config;

    static void addSinks(EnergyLedger& destination, const EnergyLedger& track) {
        destination.excitation_eV += track.excitation_eV;
        destination.ionization_eV += track.ionization_eV;
        destination.thermal_eV += track.thermal_eV;
        destination.escaped_eV += track.escaped_eV;
        destination.remaining_eV += track.remaining_eV;
    }

    static void accountCollision(EnergyLedger& ledger,
                                 const ElectronCollisionEvent& event) {
        const double lost_eV = event.energyBefore_eV - event.energyAfter_eV;
        if (!std::isfinite(lost_eV) || lost_eV < -1e-12 ||
            !std::isfinite(event.secondaryEnergy_eV) ||
            event.secondaryEnergy_eV < 0.0) {
            throw std::runtime_error("ElectronCascade: invalid collision energies");
        }
        switch (event.type) {
            case CollisionType::Elastic:
                if (std::abs(lost_eV) > 1e-12 || event.secondaryEnergy_eV != 0.0) {
                    throw std::runtime_error("ElectronCascade: elastic energy changed");
                }
                break;
            case CollisionType::Excitation:
                if (event.secondaryEnergy_eV != 0.0) {
                    throw std::runtime_error("ElectronCascade: excitation made a secondary");
                }
                ledger.excitation_eV += lost_eV;
                break;
            case CollisionType::Ionization:
                if (lost_eV + 1e-12 < event.secondaryEnergy_eV) {
                    throw std::runtime_error("ElectronCascade: ionization energy increased");
                }
                ledger.ionization_eV += lost_eV - event.secondaryEnergy_eV;
                ledger.secondaryTransfer_eV += event.secondaryEnergy_eV;
                break;
        }
    }

public:
    ElectronCascade(const ElectronTransportMC& transport,
                    ElectronCascadeConfig config)
        : transport(transport), config(config) {
        if (!std::isfinite(config.pathLength_m) || config.pathLength_m < 0.0 ||
            !std::isfinite(config.integrationStep_m) || config.integrationStep_m <= 0.0 ||
            !std::isfinite(config.thermalizationEnergy_eV) ||
            config.thermalizationEnergy_eV < 0.0 ||
            config.maxEventsPerTrack == 0 || config.maxTracks == 0) {
            throw std::invalid_argument("ElectronCascade: invalid configuration");
        }
    }

    ElectronCascadeResult run(const Vector3& start_m, const Vector3& direction,
                              double initialEnergy_eV, Random& rng) const {
        if (!start_m.isFinite() || !direction.isFinite() ||
            !std::isfinite(direction.norm()) || direction.norm() == 0.0 ||
            !std::isfinite(initialEnergy_eV) || initialEnergy_eV < 0.0) {
            throw std::invalid_argument("ElectronCascade: invalid primary electron");
        }
        ElectronCascadeResult result;
        result.ledger.initial_eV = initialEnergy_eV;
        std::deque<PendingTrack> pending;
        pending.push_back({start_m, direction / direction.norm(), initialEnergy_eV});
        double unprocessedEnergy_eV = 0.0;
        std::size_t unprocessedCount = 0;

        while (!pending.empty() && result.tracks.size() < config.maxTracks) {
            const PendingTrack input = pending.front();
            pending.pop_front();
            ElectronTrackRecord record;
            record.parentIndex = input.parentIndex;
            record.generation = input.generation;
            record.initialEnergy_eV = input.energy_eV;
            record.ledger.initial_eV = input.energy_eV;
            const std::size_t currentIndex = result.tracks.size();

            if (input.energy_eV <= config.thermalizationEnergy_eV) {
                record.thermalizedWithoutTransport = true;
                record.ledger.thermal_eV = input.energy_eV;
                record.transport = {{}, 0.0, input.energy_eV, false,
                                    input.position_m, input.direction, false};
            } else {
                record.transport = transport.run(input.position_m, input.direction,
                                                 input.energy_eV, config.pathLength_m,
                                                 config.integrationStep_m,
                                                 config.maxEventsPerTrack, rng);
                for (const ElectronCollisionEvent& event : record.transport.events) {
                    accountCollision(record.ledger, event);
                    if (event.type == CollisionType::Ionization &&
                        event.secondaryEnergy_eV > 0.0) {
                        if (currentIndex + 1 + pending.size() < config.maxTracks) {
                            pending.push_back({event.position_m, event.secondaryDirection,
                                               event.secondaryEnergy_eV, currentIndex,
                                               input.generation + 1});
                        } else {
                            unprocessedEnergy_eV += event.secondaryEnergy_eV;
                            ++unprocessedCount;
                            result.limitReached = true;
                        }
                    }
                }
                const double finalEnergy_eV = record.transport.finalEnergy_eV;
                if (record.transport.exitedAtmosphere) {
                    record.ledger.escaped_eV = finalEnergy_eV;
                } else if (finalEnergy_eV <= config.thermalizationEnergy_eV) {
                    record.ledger.thermal_eV = finalEnergy_eV;
                } else {
                    record.ledger.remaining_eV = finalEnergy_eV;
                }
                result.limitReached |= record.transport.eventLimitReached;
            }
            if (record.ledger.relativeResidual() > 1e-10) {
                throw std::runtime_error("ElectronCascade: track energy is not conserved");
            }
            result.secondaryGenerated_eV += record.ledger.secondaryTransfer_eV;
            addSinks(result.ledger, record.ledger);
            result.tracks.push_back(std::move(record));
        }

        result.unprocessedTracks = pending.size() + unprocessedCount;
        if (!pending.empty()) result.limitReached = true;
        result.ledger.remaining_eV += unprocessedEnergy_eV;
        for (const PendingTrack& input : pending) {
            result.ledger.remaining_eV += input.energy_eV;
        }
        if (result.ledger.relativeResidual() > 1e-10) {
            throw std::runtime_error("ElectronCascade: cascade energy is not conserved");
        }
        return result;
    }
};

#endif // AURORA_ELECTRONCASCADE_HPP
