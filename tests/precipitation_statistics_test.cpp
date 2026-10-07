#include "../src/experiments/PrecipitationStatistics.hpp"

#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

} // namespace

int main() {
    try {
        const auto zero = PrecipitationStatistics::wilsonInterval(0, 100);
        const auto full = PrecipitationStatistics::wilsonInterval(100, 100);
        require(zero.low == 0.0 && zero.high > 0.03 && zero.high < 0.04,
                "Wilson interval for zero events is wrong");
        require(full.high == 1.0 && full.low > 0.96 && full.low < 0.97,
                "Wilson interval for all events is wrong");
        require(std::abs(PrecipitationStatistics::theoreticalCriticalPitch_deg(0.32) -
                         69.056916384228586) < 1e-9,
                "dipole critical pitch calculation changed");

        PrecipitationStatistics::Config config{
            0.32, PrecipitationStatistics::PitchDistribution::NarrowBin,
            68.0, 0.1, 6, 12345, 1, 50};
        const auto sequential = PrecipitationStatistics::run(config);
        config.threads = 2;
        const auto parallel = PrecipitationStatistics::run(config);
        require(sequential.counts.total() == 6 &&
                parallel.counts.total() == 6 &&
                sequential.counts.precipitated == 6,
                "expected all low-pitch electrons to precipitate");
        for (std::size_t index = 0; index < 6; ++index) {
            require(sequential.initialConditions[index].pitch_rad ==
                        parallel.initialConditions[index].pitch_rad &&
                    sequential.initialConditions[index].gyroPhase_rad ==
                        parallel.initialConditions[index].gyroPhase_rad &&
                    sequential.trajectories[index].outcome ==
                        parallel.trajectories[index].outcome,
                    "parallel schedule changed a fixed-seed result");
        }
        const auto prefix = PrecipitationStatistics::count(parallel.trajectories, 3);
        require(prefix.total() == 3 && prefix.precipitated == 3,
                "outcome prefix was counted incorrectly");

        Random expectedRandom(42);
        const double u = expectedRandom.uniform01();
        config.samples = 1;
        config.seed = 42;
        config.distribution = PrecipitationStatistics::PitchDistribution::UniformAngle;
        const auto uniformAngle = PrecipitationStatistics::run(config);
        config.distribution = PrecipitationStatistics::PitchDistribution::IsotropicHemisphere;
        const auto isotropic = PrecipitationStatistics::run(config);
        require(std::abs(uniformAngle.initialConditions[0].pitch_rad -
                         u * std::numbers::pi / 2.0) < 1e-14 &&
                std::abs(isotropic.initialConditions[0].pitch_rad -
                         std::acos(1.0 - u)) < 1e-14 &&
                uniformAngle.initialConditions[0].gyroPhase_rad ==
                    isotropic.initialConditions[0].gyroPhase_rad,
                "source hypotheses sampled the wrong angular laws");
        std::cout << "Precipitation statistics passed\n";
    } catch (const std::exception& error) {
        std::cerr << "Precipitation statistics failed: " << error.what() << '\n';
        return 1;
    }
}
