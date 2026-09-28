#ifndef AURORA_SPECTRALLINECSVLOADER_HPP
#define AURORA_SPECTRALLINECSVLOADER_HPP

#include "../emission/SpectralLine.hpp"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

class SpectralLineCsvLoader {
private:
    static std::string trim(const std::string& value) {
        const auto first = value.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) return "";
        return value.substr(first, value.find_last_not_of(" \t\r\n") - first + 1);
    }

    static double number(const std::string& value) {
        std::istringstream input(value);
        double result{};
        if (!(input >> result) || !(input >> std::ws).eof() || !std::isfinite(result)) {
            throw std::invalid_argument("SpectralLineCsvLoader: invalid number " + value);
        }
        return result;
    }

public:
    static SpectralLineCatalog load(std::istream& input, const KineticsCatalog& kinetics) {
        std::string row;
        if (!std::getline(input, row) ||
            trim(row) != "id,wavelength_nm,upper_state_id,A_line_s1,source") {
            throw std::invalid_argument("SpectralLineCsvLoader: unexpected header");
        }
        std::vector<SpectralLine> lines;
        while (std::getline(input, row)) {
            if (trim(row).empty()) continue;
            std::stringstream stream(row);
            std::vector<std::string> fields;
            std::string field;
            while (std::getline(stream, field, ',')) fields.push_back(trim(field));
            if (!row.empty() && row.back() == ',') fields.push_back("");
            if (fields.size() != 5) {
                throw std::invalid_argument("SpectralLineCsvLoader: wrong number of columns");
            }
            lines.push_back({fields[0], number(fields[1]), fields[2],
                             number(fields[3]), fields[4]});
        }
        if (input.bad()) throw std::runtime_error("SpectralLineCsvLoader: read failed");
        return SpectralLineCatalog(std::move(lines), kinetics);
    }

    static SpectralLineCatalog load(const std::filesystem::path& path,
                                    const KineticsCatalog& kinetics) {
        std::ifstream input(path);
        if (!input) throw std::runtime_error("SpectralLineCsvLoader: cannot open " + path.string());
        return load(input, kinetics);
    }
};

#endif // AURORA_SPECTRALLINECSVLOADER_HPP
