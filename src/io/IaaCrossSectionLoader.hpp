#ifndef AURORA_IAACROSSSECTIONLOADER_HPP
#define AURORA_IAACROSSSECTIONLOADER_HPP

#include "../atmosphere/CrossSectionDatabase.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <istream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Reads the locally extracted IAA elmolcs v0.0.1 data. The combined file
// contains tabulated inelastic processes; integral elastic data must come
// from the species files, not from EFFECTIVE or momentum-transfer sections.
class IaaCrossSectionLoader {
private:
    struct Block {
        std::string kind;
        std::string reaction;
        std::string species;
        std::string columns;
        double threshold_eV{};
        std::vector<TabulatedPoint> points;
    };

    static std::string trim(const std::string& value) {
        const auto first = value.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) return "";
        return value.substr(first, value.find_last_not_of(" \t\r\n") - first + 1);
    }

    static bool separator(const std::string& value) {
        return value.size() >= 5 && value.find_first_not_of('-') == std::string::npos;
    }

    static bool heading(const std::string& value) {
        return value == "ELASTIC" || value == "EXCITATION" ||
               value == "IONIZATION" || value == "ROTATION" ||
               value == "ATTACHMENT" || value == "TOTAL" ||
               value == "EFFECTIVE" || value == "BREMSSTRAHLUNG";
    }

    static std::vector<Block> parse(std::istream& input, bool firstElasticOnly = false) {
        std::vector<Block> blocks;
        std::string line;
        while (std::getline(input, line)) {
            const std::string kind = trim(line);
            if (!heading(kind)) continue;
            Block block;
            block.kind = kind;
            if (!std::getline(input, block.reaction)) {
                throw std::invalid_argument("IAA: truncated process header");
            }
            block.reaction = trim(block.reaction);
            if (!std::getline(input, line)) {
                throw std::invalid_argument("IAA: missing process parameter");
            }
            if (kind != "EXCITATION" && kind != "IONIZATION" && kind != "ELASTIC") {
                continue;
            }
            std::istringstream threshold(trim(line));
            if (!(threshold >> block.threshold_eV) ||
                !std::isfinite(block.threshold_eV)) {
                throw std::invalid_argument("IAA: invalid process parameter: " + block.reaction);
            }
            bool tableStarted = false;
            bool tableEnded = false;
            while (std::getline(input, line)) {
                const std::string value = trim(line);
                if (value.rfind("SPECIES:", 0) == 0) {
                    const auto slash = value.find('/');
                    if (slash != std::string::npos) block.species = trim(value.substr(slash + 1));
                } else if (value.rfind("COLUMNS:", 0) == 0) {
                    block.columns = value;
                } else if (separator(value) && !block.columns.empty()) {
                    if (!tableStarted) tableStarted = true;
                    else { tableEnded = true; break; }
                } else if (tableStarted && !value.empty() && value[0] != '#') {
                    std::istringstream row(value);
                    TabulatedPoint point;
                    if (!(row >> point.energy_eV >> point.sigma_m2)) {
                        throw std::invalid_argument("IAA: malformed cross-section row");
                    }
                    block.points.push_back(point);
                }
            }
            if (!tableEnded || block.species.empty() || block.points.size() < 2) {
                throw std::invalid_argument("IAA: incomplete process block: " + block.reaction);
            }
            blocks.push_back(std::move(block));
            if (firstElasticOnly && kind == "ELASTIC") break;
        }
        return blocks;
    }

    static CollisionChannel makeChannel(Block block, CollisionType type,
                                        const std::string& id) {
        double scale = 1.0;
        if (block.columns.find("Energy (eV)") == std::string::npos) {
            throw std::invalid_argument("IAA: energy column is not in eV");
        }
        if (block.columns.find("Cross Section (m2)") != std::string::npos) {
            scale = 1.0;
        } else if (type == CollisionType::Elastic &&
                   block.columns.find("ICS Q0 (1e-20 m2)") != std::string::npos) {
            scale = 1e-20;
        } else {
            throw std::invalid_argument("IAA: unsupported cross-section units");
        }
        for (auto& point : block.points) point.sigma_m2 *= scale;
        double threshold_eV = 0.0;
        if (type != CollisionType::Elastic) {
            threshold_eV = block.threshold_eV;
            if (threshold_eV <= 0.0 ||
                block.points.front().energy_eV < threshold_eV * (1.0 - 1e-4)) {
                throw std::invalid_argument("IAA: invalid inelastic threshold");
            }
            if (block.points.front().energy_eV > threshold_eV) {
                block.points.insert(block.points.begin(), {threshold_eV, 0.0});
            }
        }
        return {id, block.species, type,
                type == CollisionType::Excitation ? block.reaction : "",
                threshold_eV, TabulatedFunction(std::move(block.points)),
                "IAA elmolcs v0.0.1 (Schmalzried et al., 2023), " + block.reaction};
    }

    static std::ifstream open(const std::filesystem::path& path) {
        std::ifstream input(path);
        if (!input) throw std::runtime_error("IAA: cannot open " + path.string());
        return input;
    }

public:
    static CrossSectionDatabase load(std::istream& combined, std::istream& oxygen,
                                     std::istream& oxygenMolecule,
                                     std::istream& nitrogenMolecule,
                                     double modelMaxEnergy_eV = 0.0) {
        if (!std::isfinite(modelMaxEnergy_eV) || modelMaxEnergy_eV < 0.0) {
            throw std::invalid_argument("IAA: invalid model energy limit");
        }
        std::vector<CollisionChannel> channels;
        for (Block& block : parse(combined)) {
            if (block.species != "O" && block.species != "O2" && block.species != "N2") continue;
            if (block.kind != "EXCITATION" && block.kind != "IONIZATION") continue;
            const auto type = block.kind == "EXCITATION" ? CollisionType::Excitation
                                                          : CollisionType::Ionization;
            const std::string id = "IAA:" + std::to_string(channels.size()) + ":" + block.kind;
            channels.push_back(makeChannel(std::move(block), type,
                                           id));
        }
        for (auto& [species, stream] :
             {std::pair{"O", &oxygen}, std::pair{"O2", &oxygenMolecule},
              std::pair{"N2", &nitrogenMolecule}}) {
            bool found = false;
            for (Block& block : parse(*stream, true)) {
                if (block.kind != "ELASTIC" || block.species != species) continue;
                // The first elastic block is the IAA recommended integral
                // elastic section. Later blocks are alternative datasets.
                channels.push_back(makeChannel(std::move(block), CollisionType::Elastic,
                                               std::string("IAA:") + species + ":elastic"));
                found = true;
                break;
            }
            if (!found) throw std::invalid_argument("IAA: missing integral elastic section for " +
                                                     std::string(species));
        }
        if (modelMaxEnergy_eV > 0.0) {
            for (CollisionChannel& channel : channels) {
                const double lastEnergy_eV = channel.sigma_m2.maxEnergy_eV();
                if (lastEnergy_eV >= modelMaxEnergy_eV) continue;
                // Low-energy resonances are cut off at their last knot;
                // high-energy channels use an explicitly approximate 1/E tail.
                channel.tail = lastEnergy_eV < 100.0 ? CrossSectionTail::Zero
                                                     : CrossSectionTail::InverseEnergy;
                channel.tailMaxEnergy_eV = modelMaxEnergy_eV;
                channel.source += channel.tail == CrossSectionTail::Zero
                    ? "; model: zero tail beyond measured range"
                    : "; model: inverse-energy tail beyond measured range";
            }
        }
        return CrossSectionDatabase(std::move(channels));
    }

    static CrossSectionDatabase load(const std::filesystem::path& extractedRoot,
                                     double modelMaxEnergy_eV = 0.0) {
        const auto data = extractedRoot / "elmolcs" / "Data" / "cs";
        auto combined = open(data / "Cross section.txt");
        auto oxygen = open(data / "O" / "Cross section.e-O");
        auto oxygenMolecule = open(data / "O2" / "Cross section.e-O2");
        auto nitrogenMolecule = open(data / "N2" / "Cross section.e-N2");
        return load(combined, oxygen, oxygenMolecule, nitrogenMolecule,
                    modelMaxEnergy_eV);
    }
};

#endif // AURORA_IAACROSSSECTIONLOADER_HPP
