#ifndef AURORA_LXCATCROSSSECTIONLOADER_HPP
#define AURORA_LXCATCROSSSECTIONLOADER_HPP

#include "../atmosphere/CrossSectionDatabase.hpp"

#include <cmath>
#include <fstream>
#include <istream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Imports one explicitly selected inelastic BOLSIG+/LXCat process. Selection
// is deliberate: LXCat files can contain TOTAL, momentum-transfer and partial
// cross sections that must not all be summed as independent collision events.
struct LxCatChannelSelection {
    std::string reaction;
    std::string id;
    std::string species;
    CollisionType type;
    std::string targetState;
    std::string source;
};

class LxCatCrossSectionLoader {
private:
    static std::string trim(const std::string& value) {
        const auto first = value.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) return "";
        const auto last = value.find_last_not_of(" \t\r\n");
        return value.substr(first, last - first + 1);
    }

    static bool separator(const std::string& line) {
        return line.size() >= 5 && line.find_first_not_of('-') == std::string::npos;
    }

    static std::string nextNonempty(std::istream& input) {
        std::string line;
        while (std::getline(input, line)) {
            line = trim(line);
            if (!line.empty()) return line;
        }
        throw std::invalid_argument("LXCat: truncated process block");
    }

    static double parseThreshold(const std::string& line) {
        std::istringstream stream(line);
        double threshold_eV;
        if (!(stream >> threshold_eV) || !(stream >> std::ws).eof() ||
            !std::isfinite(threshold_eV) || threshold_eV <= 0.0) {
            throw std::invalid_argument("LXCat: invalid inelastic threshold");
        }
        return threshold_eV;
    }

public:
    static CollisionChannel load(std::istream& input, const LxCatChannelSelection& selection) {
        if (selection.reaction.empty() || selection.id.empty() ||
            selection.species.empty() || selection.source.empty() ||
            (selection.type != CollisionType::Excitation &&
             selection.type != CollisionType::Ionization) ||
            (selection.type == CollisionType::Excitation && selection.targetState.empty())) {
            throw std::invalid_argument("LXCat: invalid channel selection");
        }

        const std::string wantedType = selection.type == CollisionType::Excitation
                                           ? "EXCITATION" : "IONIZATION";
        std::string line;
        std::optional<CollisionChannel> result;
        while (std::getline(input, line)) {
            if (trim(line) != wantedType) continue;
            const std::string reaction = nextNonempty(input);
            if (reaction != selection.reaction) continue;
            if (result) throw std::invalid_argument("LXCat: duplicate matching process");

            const double threshold_eV = parseThreshold(nextNonempty(input));
            bool columnsFound = false;
            bool tableStarted = false;
            std::vector<TabulatedPoint> points;
            while (std::getline(input, line)) {
                const std::string value = trim(line);
                if (value.rfind("COLUMNS:", 0) == 0) {
                    if (value.find("Energy (eV)") == std::string::npos ||
                        value.find("Cross Section (m2)") == std::string::npos) {
                        throw std::invalid_argument("LXCat: expected eV and m2 columns");
                    }
                    columnsFound = true;
                } else if (columnsFound && separator(value)) {
                    if (!tableStarted) {
                        tableStarted = true;
                    } else {
                        break;
                    }
                } else if (tableStarted && !value.empty()) {
                    std::istringstream row(value);
                    TabulatedPoint point;
                    if (!(row >> point.energy_eV >> point.sigma_m2) ||
                        !(row >> std::ws).eof()) {
                        throw std::invalid_argument("LXCat: malformed cross-section row");
                    }
                    points.push_back(point);
                }
            }
            if (!columnsFound || !tableStarted || points.empty()) {
                throw std::invalid_argument("LXCat: missing cross-section table");
            }
            if (points.front().energy_eV > threshold_eV) {
                points.insert(points.begin(), {threshold_eV, 0.0});
            }
            result.emplace(CollisionChannel{
                selection.id, selection.species, selection.type,
                selection.targetState, threshold_eV,
                TabulatedFunction(std::move(points)), selection.source});
        }
        if (!result) throw std::out_of_range("LXCat: selected process was not found");
        // Run the same metadata and threshold checks as the rest of the model.
        (void)CrossSectionDatabase({*result});
        return std::move(*result);
    }

    static CollisionChannel load(const std::string& path,
                                 const LxCatChannelSelection& selection) {
        std::ifstream input(path);
        if (!input) throw std::runtime_error("LXCat: cannot open " + path);
        return load(input, selection);
    }
};

#endif // AURORA_LXCATCROSSSECTIONLOADER_HPP
