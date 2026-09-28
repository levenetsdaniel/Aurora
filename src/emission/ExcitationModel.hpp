#ifndef AURORA_EXCITATIONMODEL_HPP
#define AURORA_EXCITATIONMODEL_HPP

#include "../atmosphere/CrossSectionDatabase.hpp"
#include "../atmosphere/ElectronCascade.hpp"
#include "../math/Vector3.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

struct ExcitedStateId {
    std::string species;
    std::string state;

    bool operator<(const ExcitedStateId& other) const {
        return std::tie(species, state) < std::tie(other.species, other.state);
    }
};

struct ExcitationGridConfig {
    Vector3 origin_m;
    Vector3 cellSize_m;
    std::size_t cellsX{};
    std::size_t cellsY{};
    std::size_t cellsZ{};
};

enum class ExcitationOutOfBounds { Reject, Ignore };

class ExcitationModel {
private:
    ExcitationGridConfig grid;
    double cellVolume_m3{};
    Vector3 upper_m;
    ExcitationOutOfBounds outOfBounds;
    std::unordered_map<std::string, ExcitedStateId> stateByChannel;
    std::map<ExcitedStateId, std::vector<double>> eventRateByState_s1;
    double ignoredOutsideRate_s1{};
    double untrackedExcitationRate_s1{};

    static std::size_t cellCount(const ExcitationGridConfig& config) {
        if (config.cellsX == 0 || config.cellsY == 0 || config.cellsZ == 0 ||
            config.cellsX > std::numeric_limits<std::size_t>::max() / config.cellsY ||
            config.cellsX * config.cellsY >
                std::numeric_limits<std::size_t>::max() / config.cellsZ) {
            throw std::invalid_argument("ExcitationModel: invalid grid dimensions");
        }
        return config.cellsX * config.cellsY * config.cellsZ;
    }

    std::size_t flatIndex(std::size_t x, std::size_t y, std::size_t z) const {
        if (x >= grid.cellsX || y >= grid.cellsY || z >= grid.cellsZ) {
            throw std::out_of_range("ExcitationModel: cell index outside grid");
        }
        return (z * grid.cellsY + y) * grid.cellsX + x;
    }

    bool eventIndex(const Vector3& position_m, std::size_t& index) const {
        if (!position_m.isFinite()) {
            throw std::invalid_argument("ExcitationModel: nonfinite event position");
        }
        if (position_m.x < grid.origin_m.x || position_m.x >= upper_m.x ||
            position_m.y < grid.origin_m.y || position_m.y >= upper_m.y ||
            position_m.z < grid.origin_m.z || position_m.z >= upper_m.z) {
            return false;
        }
        const auto x = std::min(grid.cellsX - 1, static_cast<std::size_t>(
            (position_m.x - grid.origin_m.x) / grid.cellSize_m.x));
        const auto y = std::min(grid.cellsY - 1, static_cast<std::size_t>(
            (position_m.y - grid.origin_m.y) / grid.cellSize_m.y));
        const auto z = std::min(grid.cellsZ - 1, static_cast<std::size_t>(
            (position_m.z - grid.origin_m.z) / grid.cellSize_m.z));
        index = flatIndex(x, y, z);
        return true;
    }

public:
    ExcitationModel(const CrossSectionDatabase& crossSections,
                    ExcitationGridConfig grid,
                    ExcitationOutOfBounds outOfBounds = ExcitationOutOfBounds::Reject,
                    std::vector<ExcitedStateId> selectedStates = {})
        : grid(grid), outOfBounds(outOfBounds) {
        const std::size_t count = cellCount(grid);
        if (!grid.origin_m.isFinite() || !grid.cellSize_m.isFinite() ||
            grid.cellSize_m.x <= 0.0 || grid.cellSize_m.y <= 0.0 ||
            grid.cellSize_m.z <= 0.0) {
            throw std::invalid_argument("ExcitationModel: invalid grid origin or cell size");
        }
        upper_m = {grid.origin_m.x + grid.cellSize_m.x * grid.cellsX,
                   grid.origin_m.y + grid.cellSize_m.y * grid.cellsY,
                   grid.origin_m.z + grid.cellSize_m.z * grid.cellsZ};
        cellVolume_m3 = grid.cellSize_m.x * grid.cellSize_m.y * grid.cellSize_m.z;
        if (!upper_m.isFinite() || !std::isfinite(cellVolume_m3) ||
            cellVolume_m3 <= 0.0 ||
            upper_m.x <= grid.origin_m.x || upper_m.y <= grid.origin_m.y ||
            upper_m.z <= grid.origin_m.z) {
            throw std::invalid_argument("ExcitationModel: invalid grid extent or cell volume");
        }
        if (outOfBounds != ExcitationOutOfBounds::Reject &&
            outOfBounds != ExcitationOutOfBounds::Ignore) {
            throw std::invalid_argument("ExcitationModel: invalid out-of-bounds policy");
        }
        const std::set<ExcitedStateId> selected(selectedStates.begin(), selectedStates.end());
        for (const CollisionChannel& channel : crossSections.allChannels()) {
            if (channel.type != CollisionType::Excitation) continue;
            const ExcitedStateId state{channel.species, channel.targetState};
            stateByChannel.emplace(channel.id, state);
            if (selected.empty() || selected.contains(state)) {
                eventRateByState_s1.try_emplace(state, count, 0.0);
            }
        }
        if (stateByChannel.empty() || eventRateByState_s1.empty() ||
            (!selected.empty() && selected.size() != eventRateByState_s1.size())) {
            throw std::invalid_argument("ExcitationModel: missing selected excitation channels");
        }
    }

    std::size_t accumulate(const ElectronCascadeResult& cascade,
                           double markerRate_s1) {
        if (!std::isfinite(markerRate_s1) || markerRate_s1 < 0.0) {
            throw std::invalid_argument("ExcitationModel: invalid marker rate");
        }
        if (cascade.limitReached) {
            throw std::invalid_argument("ExcitationModel: incomplete cascade");
        }
        std::map<std::pair<ExcitedStateId, std::size_t>, double> updates;
        double ignoredRate_s1 = 0.0;
        double untrackedRate_s1 = 0.0;
        std::size_t accepted = 0;
        for (const ElectronTrackRecord& track : cascade.tracks) {
            for (const ElectronCollisionEvent& event : track.transport.events) {
                if (event.type != CollisionType::Excitation) continue;
                const auto state = stateByChannel.find(event.channelId);
                if (state == stateByChannel.end()) {
                    throw std::invalid_argument("ExcitationModel: unknown excitation channel " +
                                                event.channelId);
                }
                if (!eventRateByState_s1.contains(state->second)) {
                    untrackedRate_s1 += markerRate_s1;
                    continue;
                }
                std::size_t index = 0;
                if (!eventIndex(event.position_m, index)) {
                    if (outOfBounds == ExcitationOutOfBounds::Reject) {
                        throw std::out_of_range("ExcitationModel: excitation outside grid");
                    }
                    ignoredRate_s1 += markerRate_s1;
                    continue;
                }
                updates[{state->second, index}] += markerRate_s1;
                ++accepted;
            }
        }
        if (!std::isfinite(ignoredRate_s1) ||
            !std::isfinite(ignoredOutsideRate_s1 + ignoredRate_s1) ||
            !std::isfinite(untrackedRate_s1) ||
            !std::isfinite(untrackedExcitationRate_s1 + untrackedRate_s1)) {
            throw std::overflow_error("ExcitationModel: ignored rate overflow");
        }
        for (const auto& [key, delta] : updates) {
            const double current = eventRateByState_s1.at(key.first).at(key.second);
            if (!std::isfinite(delta) || !std::isfinite(current + delta)) {
                throw std::overflow_error("ExcitationModel: event rate overflow");
            }
        }
        for (const auto& [key, delta] : updates) {
            eventRateByState_s1.at(key.first).at(key.second) += delta;
        }
        ignoredOutsideRate_s1 += ignoredRate_s1;
        untrackedExcitationRate_s1 += untrackedRate_s1;
        return accepted;
    }

    double eventRate_s1(const ExcitedStateId& state,
                        std::size_t x, std::size_t y, std::size_t z) const {
        const auto it = eventRateByState_s1.find(state);
        if (it == eventRateByState_s1.end()) {
            throw std::out_of_range("ExcitationModel: unknown excited state");
        }
        return it->second.at(flatIndex(x, y, z));
    }

    double productionRate_m3_s1(const ExcitedStateId& state,
                                 std::size_t x, std::size_t y,
                                 std::size_t z) const {
        const double rate = eventRate_s1(state, x, y, z) / cellVolume_m3;
        if (!std::isfinite(rate)) {
            throw std::overflow_error("ExcitationModel: production rate overflow");
        }
        return rate;
    }

    std::vector<ExcitedStateId> states() const {
        std::vector<ExcitedStateId> result;
        result.reserve(eventRateByState_s1.size());
        for (const auto& entry : eventRateByState_s1) {
            result.push_back(entry.first);
        }
        return result;
    }

    double ignoredEventRate_s1() const {
        return ignoredOutsideRate_s1;
    }

    double untrackedEventRate_s1() const {
        return untrackedExcitationRate_s1;
    }

    double cellVolume() const {
        return cellVolume_m3;
    }

    Vector3 cellCenter_m(std::size_t x, std::size_t y, std::size_t z) const {
        (void)flatIndex(x, y, z);
        return {grid.origin_m.x + (static_cast<double>(x) + 0.5) * grid.cellSize_m.x,
                grid.origin_m.y + (static_cast<double>(y) + 0.5) * grid.cellSize_m.y,
                grid.origin_m.z + (static_cast<double>(z) + 0.5) * grid.cellSize_m.z};
    }

    const ExcitationGridConfig& configuration() const {
        return grid;
    }
};

#endif // AURORA_EXCITATIONMODEL_HPP
