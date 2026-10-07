#ifndef AURORA_ELECTRONANIMATION_HPP
#define AURORA_ELECTRONANIMATION_HPP

#include "ElectronSceneData.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace ElectronAnimation {

struct Sample {
    Vector3 position_m;
    std::size_t completedPoints{};
    bool visible{};
    bool finished{};
};

inline void validate(const ElectronPath& path) {
    if (path.positions_m.empty() || path.positions_m.size() != path.times_s.size() ||
        !std::isfinite(path.terminalState.time_s) ||
        path.times_s.back() != path.terminalState.time_s)
        throw std::invalid_argument("ElectronAnimation: invalid recorded path");
    for (std::size_t index = 0; index < path.times_s.size(); ++index) {
        if (!std::isfinite(path.times_s[index]) || !path.positions_m[index].isFinite() ||
            (index > 0 && path.times_s[index] <= path.times_s[index - 1]))
            throw std::invalid_argument("ElectronAnimation: nonmonotone path time");
    }
}

inline Sample sample(const ElectronPath& path, double time_s) {
    if (!std::isfinite(time_s))
        throw std::invalid_argument("ElectronAnimation: invalid frame time");
    if (time_s < path.times_s.front()) return {};
    if (time_s >= path.times_s.back())
        return {path.positions_m.back(), path.positions_m.size(), true, true};
    const auto next = std::upper_bound(path.times_s.begin(), path.times_s.end(), time_s);
    const std::size_t after = static_cast<std::size_t>(next - path.times_s.begin());
    const std::size_t before = after - 1;
    const double fraction = (time_s - path.times_s[before]) /
                            (path.times_s[after] - path.times_s[before]);
    return {path.positions_m[before] * (1.0 - fraction) +
            path.positions_m[after] * fraction, after, true, false};
}

inline double lastEventTime_s(const ElectronSceneData& scene) {
    double last = 0.0;
    for (const auto& path : scene.paths) {
        validate(path);
        last = std::max(last, path.terminalState.time_s);
    }
    if (last <= 0.0) throw std::invalid_argument("ElectronAnimation: empty duration");
    return last;
}

inline double frameTime_s(std::size_t frame, std::size_t frameCount,
                          double lastTime_s) {
    if (frameCount < 2 || frame >= frameCount ||
        !std::isfinite(lastTime_s) || lastTime_s <= 0.0)
        throw std::invalid_argument("ElectronAnimation: invalid frame sequence");
    const double progress = static_cast<double>(frame) / (frameCount - 1);
    return lastTime_s * std::min(1.0, progress / 0.85);
}

inline double arrivedWeightFraction(const ElectronSceneData& scene, double time_s) {
    double total = 0.0;
    double arrived = 0.0;
    for (const auto& path : scene.paths) {
        if (path.outcome != ParticleOutcome::Precipitated) continue;
        if (!std::isfinite(path.weight_s1) || path.weight_s1 < 0.0)
            throw std::invalid_argument("ElectronAnimation: invalid particle weight");
        total += path.weight_s1;
        if (path.terminalState.time_s <= time_s) arrived += path.weight_s1;
    }
    return total > 0.0 ? arrived / total : 0.0;
}

} // namespace ElectronAnimation

#endif // AURORA_ELECTRONANIMATION_HPP
