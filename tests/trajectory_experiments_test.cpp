#include "../src/experiments/TrajectoryExperiments.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

} // namespace

int main() {
    try {
        const auto coarse = TrajectoryExperiments::uniformOrbit(50);
        const auto fine = TrajectoryExperiments::uniformOrbit(100);
        const auto finest = TrajectoryExperiments::uniformOrbit(200);
        require(coarse.positionError_m > fine.positionError_m * 3.0,
                "uniform orbit position error did not converge");
        require(fine.positionError_m > finest.positionError_m * 3.0,
                "uniform orbit position error did not converge at finer step");
        require(std::abs(fine.measuredPeriod_s - fine.exactPeriod_s) <
                std::abs(coarse.measuredPeriod_s - coarse.exactPeriod_s),
                "gyroperiod did not converge");
        require(finest.maximumRelativeEnergyError < 1e-10,
                "uniform-field energy drift is too large");

        const auto coarseEvent = TrajectoryExperiments::atmosphericCrossing(50);
        const auto fineEvent = TrajectoryExperiments::atmosphericCrossing(100);
        require(fineEvent.timeError_s < coarseEvent.timeError_s,
                "atmosphere event time did not converge");
        require(fineEvent.sphereResidual_m < 1e-3,
                "atmosphere event locator residual is too large");

        const auto mirror = TrajectoryExperiments::mirrorTrajectory(50, 200);
        require(mirror.outcome == ParticleOutcome::Mirrored,
                "control dipole particle did not mirror");
        require(mirror.steps > 0 && mirror.trace.size() > 2,
                "mirror trace is empty");
        require(mirror.mirrorMagneticMagnitude_T > mirror.initialMagneticMagnitude_T,
                "mirror did not occur in a stronger field");
        require(std::abs(mirror.mirrorMagneticMagnitude_T /
                         mirror.predictedMirrorMagneticMagnitude_T - 1.0) < 0.05,
                "mirror field differs too much from adiabatic prediction");
        require(std::abs(mirror.relativeEnergyError) < 1e-8,
                "dipole energy drift is too large");
        std::cout << "Trajectory experiments passed\n";
    } catch (const std::exception& error) {
        std::cerr << "Trajectory experiments failed: " << error.what() << '\n';
        return 1;
    }
}
