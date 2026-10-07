#ifndef AURORA_EMISSIONCSVREPLAY_HPP
#define AURORA_EMISSIONCSVREPLAY_HPP

#include "../emission/VolumeEmissionGrid.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace EmissionCsvReplay {

inline std::vector<std::string> fields(const std::string& line) {
    std::vector<std::string> result;
    std::istringstream input(line);
    std::string field;
    while (std::getline(input, field, ',')) result.push_back(field);
    return result;
}

inline VolumeEmissionGrid load(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open emission CSV: " + path.string());
    std::string line;
    const auto metadata = [&](const std::string& prefix) {
        if (!std::getline(input, line) || line.rfind(prefix, 0) != 0)
            throw std::invalid_argument("Emission CSV missing " + prefix);
        return fields(line.substr(prefix.size()));
    };
    const auto originFields = metadata("# grid_origin_m=");
    const auto sizeFields = metadata("# cell_size_m=");
    const auto shapeFields = metadata("# shape_xyz=");
    if (originFields.size() != 3 || sizeFields.size() != 3 ||
        shapeFields.size() != 3 || !std::getline(input, line) ||
        line != "x_m,y_m,z_m,j_O_5577_ph_m3_s1,j_O_6300_ph_m3_s1")
        throw std::invalid_argument("Unsupported emission CSV layout");
    const Vector3 origin(std::stod(originFields[0]), std::stod(originFields[1]),
                         std::stod(originFields[2]));
    const Vector3 size(std::stod(sizeFields[0]), std::stod(sizeFields[1]),
                       std::stod(sizeFields[2]));
    const ExcitationGridConfig config{origin, size, std::stoull(shapeFields[0]),
                                      std::stoull(shapeFields[1]),
                                      std::stoull(shapeFields[2])};
    VolumeEmissionGrid grid(config, {{"O_5577", 557.7, "replay-green", 1.0, "CSV replay"},
                                     {"O_6300", 630.0, "replay-red", 1.0, "CSV replay"}});
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto values = fields(line);
        if (values.size() != 5)
            throw std::invalid_argument("Malformed emission CSV row");
        const std::array<double, 3> coordinate{
            std::stod(values[0]), std::stod(values[1]), std::stod(values[2])};
        const std::array<double, 3> lower{origin.x, origin.y, origin.z};
        const std::array<double, 3> step{size.x, size.y, size.z};
        const std::array<std::size_t, 3> shape{config.cellsX, config.cellsY,
                                               config.cellsZ};
        std::array<std::size_t, 3> cell{};
        for (std::size_t axis = 0; axis < 3; ++axis) {
            const double raw = (coordinate[axis] - lower[axis]) / step[axis] - 0.5;
            if (!std::isfinite(raw) || raw < 0.0 ||
                raw >= static_cast<double>(shape[axis]))
                throw std::invalid_argument("Emission CSV point outside grid");
            const long long rounded = std::llround(raw);
            if (rounded < 0 || rounded >= static_cast<long long>(shape[axis]) ||
                std::abs(raw - rounded) > 1e-6)
                throw std::invalid_argument("Emission CSV point outside grid");
            cell[axis] = static_cast<std::size_t>(rounded);
        }
        grid.setEmissivity("O_5577", cell[0], cell[1], cell[2], std::stod(values[3]));
        grid.setEmissivity("O_6300", cell[0], cell[1], cell[2], std::stod(values[4]));
    }
    if (!input.eof()) throw std::runtime_error("Emission CSV read failed");
    return grid;
}

inline void validateCompanionManifest(const std::filesystem::path& csv,
                                      std::uint64_t seed, std::size_t markers,
                                      double sourceRate_s1,
                                      double minPitch_deg, double maxPitch_deg) {
    std::string stem = csv.stem().string();
    constexpr const char* suffix = "_emission";
    if (stem.size() <= std::char_traits<char>::length(suffix) ||
        stem.substr(stem.size() - std::char_traits<char>::length(suffix)) != suffix)
        throw std::invalid_argument("Replay CSV name must end in _emission.csv");
    stem.resize(stem.size() - std::char_traits<char>::length(suffix));
    const auto manifestPath = csv.parent_path() / (stem + "_manifest.txt");
    std::ifstream manifest(manifestPath);
    if (!manifest) throw std::runtime_error("Replay requires companion manifest: " +
                                            manifestPath.string());
    bool seedMatched = false;
    bool markersMatched = false;
    bool rateMatched = false;
    double storedMinPitch_deg = 0.0;
    double storedMaxPitch_deg = 90.0;
    std::string line;
    while (std::getline(manifest, line)) {
        if (line.rfind("seed=", 0) == 0)
            seedMatched = std::stoull(line.substr(5)) == seed;
        else if (line.rfind("visual_markers=", 0) == 0)
            markersMatched = std::stoull(line.substr(15)) == markers;
        else if (line.rfind("source_rate_electrons_s1=", 0) == 0) {
            const double value = std::stod(line.substr(25));
            rateMatched = std::abs(value / sourceRate_s1 - 1.0) < 1e-9;
        } else if (line.rfind("pitch_min_deg=", 0) == 0)
            storedMinPitch_deg = std::stod(line.substr(14));
        else if (line.rfind("pitch_max_deg=", 0) == 0)
            storedMaxPitch_deg = std::stod(line.substr(14));
    }
    if (!seedMatched || !markersMatched || !rateMatched ||
        std::abs(storedMinPitch_deg - minPitch_deg) > 1e-9 ||
        std::abs(storedMaxPitch_deg - maxPitch_deg) > 1e-9)
        throw std::invalid_argument("Replay source configuration differs from manifest");
}

} // namespace EmissionCsvReplay

#endif // AURORA_EMISSIONCSVREPLAY_HPP
