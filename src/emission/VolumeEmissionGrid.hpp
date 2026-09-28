#ifndef AURORA_VOLUMEEMISSIONGRID_HPP
#define AURORA_VOLUMEEMISSIONGRID_HPP

#include "SpectralLine.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Cell values are photon production rates in photons m^-3 s^-1.
class VolumeEmissionGrid {
private:
    ExcitationGridConfig grid;
    Vector3 upper_m;
    std::vector<SpectralLine> spectralLines;
    std::map<std::string, std::vector<double>> emissivityByLine;

    std::size_t index(std::size_t x, std::size_t y, std::size_t z) const {
        if (x >= grid.cellsX || y >= grid.cellsY || z >= grid.cellsZ) {
            throw std::out_of_range("VolumeEmissionGrid: cell outside grid");
        }
        return (z * grid.cellsY + y) * grid.cellsX + x;
    }

public:
    VolumeEmissionGrid(ExcitationGridConfig config, std::vector<SpectralLine> lines)
        : grid(config), spectralLines(std::move(lines)) {
        if (!grid.origin_m.isFinite() || !grid.cellSize_m.isFinite() ||
            grid.cellSize_m.x <= 0.0 || grid.cellSize_m.y <= 0.0 ||
            grid.cellSize_m.z <= 0.0 || grid.cellsX == 0 || grid.cellsY == 0 ||
            grid.cellsZ == 0 || spectralLines.empty() ||
            grid.cellsX > std::numeric_limits<std::size_t>::max() / grid.cellsY ||
            grid.cellsX * grid.cellsY >
                std::numeric_limits<std::size_t>::max() / grid.cellsZ) {
            throw std::invalid_argument("VolumeEmissionGrid: invalid geometry or lines");
        }
        upper_m = {grid.origin_m.x + grid.cellSize_m.x * grid.cellsX,
                   grid.origin_m.y + grid.cellSize_m.y * grid.cellsY,
                   grid.origin_m.z + grid.cellSize_m.z * grid.cellsZ};
        if (!upper_m.isFinite() || upper_m.x <= grid.origin_m.x ||
            upper_m.y <= grid.origin_m.y || upper_m.z <= grid.origin_m.z) {
            throw std::invalid_argument("VolumeEmissionGrid: invalid extent");
        }
        const std::size_t count = grid.cellsX * grid.cellsY * grid.cellsZ;
        for (const SpectralLine& line : spectralLines) {
            if (line.id.empty() || line.upperStateId.empty() || line.source.empty() ||
                !std::isfinite(line.wavelength_nm) || line.wavelength_nm <= 0.0 ||
                !std::isfinite(line.transitionRate_s1) || line.transitionRate_s1 <= 0.0 ||
                !emissivityByLine.emplace(line.id, std::vector<double>(count, 0.0)).second) {
                throw std::invalid_argument("VolumeEmissionGrid: invalid or duplicate line");
            }
        }
    }

    static VolumeEmissionGrid fromKinetics(const ExcitationModel& excitation,
                                            const ExcitedStateKinetics& kinetics,
                                            const KineticsCatalog& states,
                                            const SpectralLineCatalog& lines) {
        VolumeEmissionGrid result(excitation.configuration(), lines.lines());
        const std::size_t count = result.grid.cellsX * result.grid.cellsY * result.grid.cellsZ;
        for (const ExcitedStateId& eventState : excitation.states()) {
            bool produced = false;
            for (std::size_t z = 0; z < result.grid.cellsZ && !produced; ++z) {
                for (std::size_t y = 0; y < result.grid.cellsY && !produced; ++y) {
                    for (std::size_t x = 0; x < result.grid.cellsX; ++x) {
                        if (excitation.eventRate_s1(eventState, x, y, z) > 0.0) {
                            produced = true;
                            break;
                        }
                    }
                }
            }
            if (!produced) continue;
            const RadiativeState& state = states.resolve(eventState);
            std::vector<const SpectralLine*> matchingLines;
            for (const SpectralLine& line : lines.lines()) {
                if (line.upperStateId == state.id) matchingLines.push_back(&line);
            }
            if (matchingLines.empty()) {
                throw std::out_of_range("VolumeEmissionGrid: no line for produced state " +
                                        state.id);
            }
            const std::vector<StatePopulation> populations =
                kinetics.solveGrid(excitation, eventState);
            if (populations.size() != count) {
                throw std::logic_error("VolumeEmissionGrid: population size mismatch");
            }
            for (const SpectralLine* line : matchingLines) {
                std::vector<double>& values = result.emissivityByLine.at(line->id);
                for (std::size_t cell = 0; cell < count; ++cell) {
                    const double contribution =
                        line->transitionRate_s1 * populations[cell].density_m3;
                    const double updated = values[cell] + contribution;
                    if (!std::isfinite(updated) || updated < 0.0) {
                        throw std::overflow_error("VolumeEmissionGrid: emissivity overflow");
                    }
                    values[cell] = updated;
                }
            }
        }
        return result;
    }

    void setEmissivity(const std::string& lineId,
                       std::size_t x, std::size_t y, std::size_t z,
                       double value_ph_m3_s1) {
        if (!std::isfinite(value_ph_m3_s1) || value_ph_m3_s1 < 0.0) {
            throw std::invalid_argument("VolumeEmissionGrid: invalid emissivity");
        }
        emissivityByLine.at(lineId).at(index(x, y, z)) = value_ph_m3_s1;
    }

    double emissivity(const std::string& lineId,
                      std::size_t x, std::size_t y, std::size_t z) const {
        return emissivityByLine.at(lineId).at(index(x, y, z));
    }

    const ExcitationGridConfig& configuration() const { return grid; }
    const Vector3& upperCorner_m() const { return upper_m; }
    const std::vector<SpectralLine>& lines() const { return spectralLines; }
};

#endif // AURORA_VOLUMEEMISSIONGRID_HPP
