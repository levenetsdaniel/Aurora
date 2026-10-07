#include "../src/visualization/ElectronAnimation.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

int main() {
    try {
        ElectronPath precipitating;
        precipitating.weight_s1 = 3.0;
        precipitating.outcome = ParticleOutcome::Precipitated;
        precipitating.positions_m = {{1, 0, 0}, {3, 0, 0}, {5, 0, 0}};
        precipitating.times_s = {0.0, 1.0, 2.0};
        precipitating.terminalState.time_s = 2.0;
        ElectronPath mirrored;
        mirrored.weight_s1 = 1.0;
        mirrored.outcome = ParticleOutcome::Mirrored;
        mirrored.positions_m = {{1, 1, 0}, {1, 3, 0}};
        mirrored.times_s = {0.0, 4.0};
        mirrored.terminalState.time_s = 4.0;
        const ElectronSceneData scene{150000.0, {precipitating, mirrored}, {}};
        const auto moving = ElectronAnimation::sample(scene.paths[0], 0.5);
        const auto finished = ElectronAnimation::sample(scene.paths[0], 3.0);
        if (!moving.visible || moving.finished || moving.completedPoints != 1 ||
            std::abs(moving.position_m.x - 2.0) > 1e-12 ||
            !finished.finished || finished.position_m.x != 5.0 ||
            ElectronAnimation::lastEventTime_s(scene) != 4.0 ||
            ElectronAnimation::arrivedWeightFraction(scene, 1.99) != 0.0 ||
            ElectronAnimation::arrivedWeightFraction(scene, 2.0) != 1.0 ||
            ElectronAnimation::frameTime_s(0, 21, 4.0) != 0.0 ||
            ElectronAnimation::frameTime_s(20, 21, 4.0) != 4.0)
            throw std::runtime_error("Time interpolation or emission onset is wrong");
        precipitating.times_s[1] = 0.0;
        bool rejected = false;
        try { ElectronAnimation::validate(precipitating); }
        catch (const std::invalid_argument&) { rejected = true; }
        if (!rejected) throw std::runtime_error("Nonmonotone time was accepted");
        std::cout << "Electron animation checks passed\n";
    } catch (const std::exception& error) {
        std::cerr << "Electron animation failed: " << error.what() << '\n';
        return 1;
    }
}
