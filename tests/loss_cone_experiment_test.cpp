#include "../src/experiments/LossConeExperiments.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

int main() {
    try {
        const auto result = LossConeExperiments::run(0.32, 150000.0, 50, 0.5);
        if (result.L <= 1.0 ||
            result.precipitatingPitch_deg >= result.mirroringPitch_deg ||
            result.mirroringPitch_deg - result.precipitatingPitch_deg > 0.5 ||
            result.boundaryField_T <= result.initialField_T ||
            std::abs(result.numericalPitch_deg - result.theoreticalPitch_deg) > 2.0) {
            throw std::runtime_error("loss-cone boundary does not match dipole estimate");
        }
        try {
            (void)LossConeExperiments::run(0.32, 450000.0, 50);
            throw std::runtime_error("invalid boundary was accepted");
        } catch (const std::invalid_argument&) {
        }
        std::cout << "Loss-cone experiment passed\n";
    } catch (const std::exception& error) {
        std::cerr << "Loss-cone experiment failed: " << error.what() << '\n';
        return 1;
    }
}
