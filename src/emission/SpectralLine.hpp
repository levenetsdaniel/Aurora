#ifndef AURORA_SPECTRALLINE_HPP
#define AURORA_SPECTRALLINE_HPP

#include "ExcitedStateKinetics.hpp"

#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

struct SpectralLine {
    std::string id;
    double wavelength_nm{};
    std::string upperStateId;
    double transitionRate_s1{};
    std::string source;
};

class SpectralLineCatalog {
private:
    std::vector<SpectralLine> spectralLines;

public:
    SpectralLineCatalog(std::vector<SpectralLine> lines,
                        const KineticsCatalog& kinetics) : spectralLines(std::move(lines)) {
        if (spectralLines.empty()) {
            throw std::invalid_argument("SpectralLineCatalog: no lines");
        }
        std::map<std::string, bool> seenIds;
        std::map<std::string, double> ratesByState;
        for (const SpectralLine& line : spectralLines) {
            if (line.id.empty() || line.upperStateId.empty() || line.source.empty() ||
                !std::isfinite(line.wavelength_nm) || line.wavelength_nm <= 0.0 ||
                !std::isfinite(line.transitionRate_s1) || line.transitionRate_s1 <= 0.0 ||
                !seenIds.emplace(line.id, true).second) {
                throw std::invalid_argument("SpectralLineCatalog: invalid or duplicate line");
            }
            ratesByState[line.upperStateId] += line.transitionRate_s1;
        }
        // Each A_lambda is an absolute transition rate, never a branching ratio.
        for (const auto& [stateId, sum_s1] : ratesByState) {
            const RadiativeState& state = kinetics.stateById(stateId);
            if (!std::isfinite(sum_s1) ||
                sum_s1 > state.totalDecay_s1 * (1.0 + 1e-10)) {
                throw std::invalid_argument("SpectralLineCatalog: line rates exceed A_total for " +
                                            stateId);
            }
        }
    }

    const std::vector<SpectralLine>& lines() const { return spectralLines; }
};

#endif // AURORA_SPECTRALLINE_HPP
