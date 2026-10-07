#include "../src/visualization/EmissionCsvReplay.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

int main() {
    try {
        const auto path = std::filesystem::path(__FILE__).parent_path() /
                          "fixtures/replay_emission.csv";
        EmissionCsvReplay::validateCompanionManifest(path, 42, 2, 1e20, 0, 90);
        const auto grid = EmissionCsvReplay::load(path);
        if (grid.configuration().cellsX != 2 ||
            grid.emissivity("O_5577", 0, 0, 0) != 10.0 ||
            grid.emissivity("O_6300", 1, 0, 0) != 40.0)
            throw std::runtime_error("Emission replay changed voxel values");
        bool rejected = false;
        try { EmissionCsvReplay::validateCompanionManifest(path, 43, 2, 1e20, 0, 90); }
        catch (const std::invalid_argument&) { rejected = true; }
        if (!rejected) throw std::runtime_error("Mismatched seed was accepted");
        rejected = false;
        try { EmissionCsvReplay::validateCompanionManifest(path, 42, 2, 1e20, 10, 90); }
        catch (const std::invalid_argument&) { rejected = true; }
        if (!rejected) throw std::runtime_error("Mismatched pitch was accepted");
        std::cout << "Emission CSV replay checks passed\n";
    } catch (const std::exception& error) {
        std::cerr << "Emission CSV replay failed: " << error.what() << '\n';
        return 1;
    }
}
