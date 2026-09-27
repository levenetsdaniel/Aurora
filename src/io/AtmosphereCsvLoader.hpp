#ifndef AURORA_ATMOSPHERECSVLOADER_HPP
#define AURORA_ATMOSPHERECSVLOADER_HPP

#include "../atmosphere/AtmosphereModel.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <istream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

class AtmosphereCsvLoader {
private:
    static std::string trim(const std::string& value) {
        const std::size_t first = value.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) return {};
        const std::size_t last = value.find_last_not_of(" \t\r\n");
        return value.substr(first, last - first + 1);
    }

    static AtmosphereSample parseRow(const std::string& line, std::size_t lineNumber) {
        if (std::count(line.begin(), line.end(), ',') != 4) {
            throw std::runtime_error("Atmosphere CSV: expected five columns at line " +
                                     std::to_string(lineNumber));
        }
        std::istringstream stream(line);
        std::array<double, 5> values{};
        std::string token;
        for (double& value : values) {
            if (!std::getline(stream, token, ',')) {
                throw std::runtime_error("Atmosphere CSV: missing column at line " +
                                         std::to_string(lineNumber));
            }
            token = trim(token);
            std::size_t consumed = 0;
            try {
                value = std::stod(token, &consumed);
            } catch (const std::exception&) {
                throw std::runtime_error("Atmosphere CSV: invalid number at line " +
                                         std::to_string(lineNumber));
            }
            if (consumed != token.size()) {
                throw std::runtime_error("Atmosphere CSV: invalid number at line " +
                                         std::to_string(lineNumber));
            }
        }
        if (std::getline(stream, token, ',')) {
            throw std::runtime_error("Atmosphere CSV: extra column at line " +
                                     std::to_string(lineNumber));
        }
        return {values[0], values[1], values[2], values[3], values[4]};
    }

public:
    static AtmosphereModel parse(std::istream& input) {
        constexpr const char* expectedHeader =
            "altitude_m,nO_m3,nO2_m3,nN2_m3,temperature_K";
        std::vector<AtmosphereSample> samples;
        std::string sourceMetadata;
        bool sourceFound = false;
        bool headerFound = false;
        std::string line;
        std::size_t lineNumber = 0;

        while (std::getline(input, line)) {
            ++lineNumber;
            line = trim(line);
            if (line.empty()) continue;
            if (line.front() == '#') {
                if (headerFound) {
                    throw std::runtime_error("Atmosphere CSV: metadata must precede the header");
                }
                const std::string metadata = trim(line.substr(1));
                if (metadata.rfind("source=", 0) == 0 && metadata.size() > 7) {
                    sourceFound = true;
                }
                sourceMetadata += metadata + '\n';
                continue;
            }
            if (!headerFound) {
                if (line != expectedHeader) {
                    throw std::runtime_error("Atmosphere CSV: expected SI-unit column header");
                }
                headerFound = true;
                continue;
            }
            samples.push_back(parseRow(line, lineNumber));
        }

        if (input.bad()) {
            throw std::runtime_error("Atmosphere CSV: error while reading input");
        }
        if (!headerFound || !sourceFound) {
            throw std::runtime_error("Atmosphere CSV: missing header or source metadata");
        }
        return AtmosphereModel(std::move(samples), std::move(sourceMetadata));
    }

    static AtmosphereModel load(const std::filesystem::path& path) {
        std::ifstream input(path);
        if (!input) {
            throw std::runtime_error("Atmosphere CSV: cannot open " + path.string());
        }
        return parse(input);
    }
};

#endif // AURORA_ATMOSPHERECSVLOADER_HPP
