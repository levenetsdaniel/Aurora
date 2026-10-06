#include "LossConeExperiments.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main(int argc, char* argv[]) {
    try {
        if (argc != 2) {
            std::cerr << "Usage: aurora_loss_cone_experiment OUTPUT_DIRECTORY\n";
            return 2;
        }
        const std::filesystem::path directory(argv[1]);
        std::filesystem::create_directories(directory);
        std::ofstream output(directory / "loss_cone.csv");
        if (!output) throw std::runtime_error("Cannot write loss_cone.csv");
        output << std::setprecision(17);
        output << "source_colatitude_deg,L,atmosphere_height_m,steps_per_gyroperiod,"
                  "dt_s,precipitating_pitch_deg,mirroring_pitch_deg,"
                  "numerical_pitch_deg,theoretical_pitch_deg,B_initial_T,"
                  "B_boundary_T,iterations\n";
        const auto write = [&](double colatitude, double atmosphereHeight_m,
                               int stepsPerGyroperiod) {
            std::cout << "Scanning colatitude " << colatitude * 180.0 /
                         std::numbers::pi << " deg, boundary " << atmosphereHeight_m /
                         1000.0 << " km, T/dt=" << stepsPerGyroperiod << std::endl;
            const auto result = LossConeExperiments::run(
                colatitude, atmosphereHeight_m, stepsPerGyroperiod);
            output << result.sourceColatitude_deg << ',' << result.L << ','
                   << result.atmosphereHeight_m << ',' << result.stepsPerGyroperiod << ','
                   << result.dt_s << ',' << result.precipitatingPitch_deg << ','
                   << result.mirroringPitch_deg << ',' << result.numericalPitch_deg << ','
                   << result.theoreticalPitch_deg << ',' << result.initialField_T << ','
                   << result.boundaryField_T << ',' << result.iterations << '\n';
            std::cout << "  numerical " << result.numericalPitch_deg
                      << " deg; theory " << result.theoreticalPitch_deg << " deg\n";
        };
        write(0.25, 150000.0, 50);
        write(0.32, 150000.0, 50);
        write(0.40, 150000.0, 50);
        write(0.32, 150000.0, 100);
        write(0.32, 100000.0, 50);
        write(0.32, 200000.0, 50);
        if (!output) throw std::runtime_error("Failed writing loss_cone.csv");
        std::cout << "Saved " << directory / "loss_cone.csv" << '\n';
    } catch (const std::exception& error) {
        std::cerr << "Loss-cone experiment failed: " << error.what() << '\n';
        return 1;
    }
}
