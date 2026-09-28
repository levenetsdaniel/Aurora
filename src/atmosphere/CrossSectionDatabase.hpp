#ifndef AURORA_CROSSSECTIONDATABASE_HPP
#define AURORA_CROSSSECTIONDATABASE_HPP

#include "../numerics/Interpolation.hpp"

#include <cmath>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

enum class CollisionType { Elastic, Excitation, Ionization };
enum class CrossSectionTail { Reject, Zero, InverseEnergy };

struct CollisionChannel {
    std::string id;
    std::string species;
    CollisionType type;
    std::string targetState;
    double threshold_eV{};
    TabulatedFunction sigma_m2;
    std::string source;
    CrossSectionTail tail{CrossSectionTail::Reject};
    double tailMaxEnergy_eV{};
};

// Inelastic tables begin at their threshold with sigma=0. Elastic tables may
// begin above 0 eV; queries below that first knot are rejected as out of range.
class CrossSectionDatabase {
private:
    std::vector<CollisionChannel> channels;

    static void validateEnergy(double energy_eV) {
        if (!std::isfinite(energy_eV) || energy_eV < 0.0) {
            throw std::invalid_argument("CrossSectionDatabase: energy must be finite and nonnegative");
        }
    }

public:
    explicit CrossSectionDatabase(std::vector<CollisionChannel> channels)
        : channels(std::move(channels)) {
        if (this->channels.empty()) {
            throw std::invalid_argument("CrossSectionDatabase: at least one channel is required");
        }

        std::unordered_set<std::string> ids;
        for (const CollisionChannel& channel : this->channels) {
            if (channel.id.empty() || !ids.insert(channel.id).second ||
                channel.species.empty() || channel.source.empty() ||
                !std::isfinite(channel.threshold_eV) || channel.threshold_eV < 0.0) {
                throw std::invalid_argument("CrossSectionDatabase: invalid channel metadata or threshold");
            }
            switch (channel.type) {
                case CollisionType::Elastic:
                    if (channel.threshold_eV != 0.0) {
                        throw std::invalid_argument("CrossSectionDatabase: elastic threshold must be zero");
                    }
                    break;
                case CollisionType::Excitation:
                    if (channel.threshold_eV <= 0.0 || channel.targetState.empty() ||
                        channel.sigma_m2.minEnergy_eV() != channel.threshold_eV) {
                        throw std::invalid_argument("CrossSectionDatabase: excitation needs a state and positive threshold");
                    }
                    break;
                case CollisionType::Ionization:
                    if (channel.threshold_eV <= 0.0 ||
                        channel.sigma_m2.minEnergy_eV() != channel.threshold_eV) {
                        throw std::invalid_argument("CrossSectionDatabase: ionization needs a positive threshold");
                    }
                    break;
                default:
                    throw std::invalid_argument("CrossSectionDatabase: unknown collision type");
            }
            if (channel.type != CollisionType::Elastic &&
                channel.sigma_m2.firstValue_m2() != 0.0) {
                throw std::invalid_argument("CrossSectionDatabase: inelastic cross section must vanish at threshold");
            }
            if (channel.tail == CrossSectionTail::Reject) {
                if (channel.tailMaxEnergy_eV != 0.0) {
                    throw std::invalid_argument("CrossSectionDatabase: reject tail cannot have an energy limit");
                }
            } else if (!std::isfinite(channel.tailMaxEnergy_eV) ||
                       channel.tailMaxEnergy_eV <= channel.sigma_m2.maxEnergy_eV()) {
                throw std::invalid_argument("CrossSectionDatabase: invalid tail energy limit");
            }
        }
    }

    double sigma(const CollisionChannel& channel, double energy_eV) const {
        validateEnergy(energy_eV);
        if (energy_eV < channel.threshold_eV) return 0.0;
        if (energy_eV > channel.sigma_m2.maxEnergy_eV()) {
            if (channel.tail == CrossSectionTail::Reject ||
                energy_eV > channel.tailMaxEnergy_eV) {
                throw std::out_of_range("CrossSectionDatabase: channel " + channel.id +
                                        " does not cover requested energy");
            }
            if (channel.tail == CrossSectionTail::Zero) return 0.0;
            return channel.sigma_m2.lastValue_m2() *
                channel.sigma_m2.maxEnergy_eV() / energy_eV;
        }
        return channel.sigma_m2.evaluate(energy_eV);
    }

    double totalSigma(const std::string& species, double energy_eV) const {
        validateEnergy(energy_eV);
        if (species.empty()) {
            throw std::invalid_argument("CrossSectionDatabase: species must not be empty");
        }

        bool found = false;
        double total_m2 = 0.0;
        for (const CollisionChannel& channel : channels) {
            if (channel.species == species) {
                found = true;
                total_m2 += sigma(channel, energy_eV);
            }
        }
        if (!found) {
            throw std::out_of_range("CrossSectionDatabase: species has no channels");
        }
        if (!std::isfinite(total_m2)) {
            throw std::overflow_error("CrossSectionDatabase: total cross section overflow");
        }
        return total_m2;
    }

    // Checks tabulated energy coverage only, not whether every physical
    // collision process is represented by a channel.
    void requireEnergyCoverage(const std::string& species,
                               double minEnergy_eV, double maxEnergy_eV) const {
        validateEnergy(minEnergy_eV);
        validateEnergy(maxEnergy_eV);
        if (species.empty() || maxEnergy_eV < minEnergy_eV) {
            throw std::invalid_argument("CrossSectionDatabase: invalid coverage request");
        }
        bool found = false;
        for (const CollisionChannel& channel : channels) {
            if (channel.species != species) continue;
            found = true;
            if (maxEnergy_eV > (channel.tail == CrossSectionTail::Reject
                                    ? channel.sigma_m2.maxEnergy_eV()
                                    : channel.tailMaxEnergy_eV) ||
                (channel.type == CollisionType::Elastic &&
                 minEnergy_eV < channel.sigma_m2.minEnergy_eV())) {
                throw std::out_of_range(
                    "CrossSectionDatabase: channel " + channel.id +
                    " does not cover requested energies");
            }
        }
        if (!found) {
            throw std::out_of_range("CrossSectionDatabase: species has no channels");
        }
    }

    std::vector<const CollisionChannel*> activeChannels(double energy_eV) const {
        validateEnergy(energy_eV);
        std::vector<const CollisionChannel*> active;
        for (const CollisionChannel& channel : channels) {
            if (sigma(channel, energy_eV) > 0.0) active.push_back(&channel);
        }
        return active;
    }

    const CollisionChannel& channelById(std::string_view id) const {
        for (const CollisionChannel& channel : channels) {
            if (channel.id == id) return channel;
        }
        throw std::out_of_range("CrossSectionDatabase: channel id not found");
    }

    const std::vector<CollisionChannel>& allChannels() const {
        return channels;
    }
};

#endif // AURORA_CROSSSECTIONDATABASE_HPP
