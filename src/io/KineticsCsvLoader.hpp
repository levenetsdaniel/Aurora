#ifndef AURORA_KINETICSCSVLOADER_HPP
#define AURORA_KINETICSCSVLOADER_HPP

#include "../emission/ExcitedStateKinetics.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <istream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Strict, unquoted CSV format for the versioned kinetics tables in data/kinetics.
// Labels and URLs in these tables must not contain commas.
class KineticsCsvLoader {
private:
    static std::string trim(const std::string& value) {
        const auto first = value.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) return "";
        return value.substr(first, value.find_last_not_of(" \t\r\n") - first + 1);
    }

    static std::vector<std::string> fields(const std::string& line, std::size_t count) {
        std::vector<std::string> result;
        std::stringstream stream(line);
        std::string field;
        while (std::getline(stream, field, ',')) result.push_back(trim(field));
        if (!line.empty() && line.back() == ',') result.push_back("");
        if (result.size() != count) {
            throw std::invalid_argument("KineticsCsvLoader: wrong number of columns");
        }
        return result;
    }

    static double number(const std::string& value) {
        std::istringstream stream(value);
        double result;
        if (!(stream >> result) || !(stream >> std::ws).eof() ||
            !std::isfinite(result)) {
            throw std::invalid_argument("KineticsCsvLoader: invalid number " + value);
        }
        return result;
    }

    static void header(std::istream& input, const std::string& expected) {
        std::string line;
        if (!std::getline(input, line) || trim(line) != expected) {
            throw std::invalid_argument("KineticsCsvLoader: unexpected CSV header");
        }
    }

    static std::ifstream open(const std::filesystem::path& path) {
        std::ifstream input(path);
        if (!input) throw std::runtime_error("KineticsCsvLoader: cannot open " + path.string());
        return input;
    }

public:
    static KineticsCatalog load(std::istream& statesInput,
                                std::istream& quenchingInput,
                                std::istream& aliasesInput) {
        header(statesInput, "id,species,A_total_s1,source");
        header(quenchingInput,
               "id,quencher,k0_m3_s,T_ref_K,beta,activation_K,source");
        header(aliasesInput, "species,event_state,id");
        std::vector<RadiativeState> states;
        std::vector<QuenchingTerm> terms;
        std::vector<StateAlias> aliases;
        std::string line;
        while (std::getline(statesInput, line)) {
            if (trim(line).empty()) continue;
            const auto row = fields(line, 4);
            states.push_back({row[0], row[1], number(row[2]), row[3]});
        }
        while (std::getline(quenchingInput, line)) {
            if (trim(line).empty()) continue;
            const auto row = fields(line, 7);
            terms.push_back({row[0], row[1], number(row[2]), number(row[3]),
                             number(row[4]), number(row[5]), row[6]});
        }
        while (std::getline(aliasesInput, line)) {
            if (trim(line).empty()) continue;
            const auto row = fields(line, 3);
            aliases.push_back({{row[0], row[1]}, row[2]});
        }
        if (statesInput.bad() || quenchingInput.bad() || aliasesInput.bad()) {
            throw std::runtime_error("KineticsCsvLoader: failed to read CSV");
        }
        return KineticsCatalog(std::move(states), std::move(terms),
                               std::move(aliases));
    }

    static KineticsCatalog load(const std::filesystem::path& directory) {
        auto states = open(directory / "states.csv");
        auto quenching = open(directory / "quenching.csv");
        auto aliases = open(directory / "aliases.csv");
        return load(states, quenching, aliases);
    }
};

#endif // AURORA_KINETICSCSVLOADER_HPP
