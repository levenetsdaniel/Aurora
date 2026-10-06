#ifndef AURORA_PRECIPITATIONSTATISTICS_HPP
#define AURORA_PRECIPITATIONSTATISTICS_HPP

#include "TrajectoryExperiments.hpp"
#include "../ensemble/InitialCondition.hpp"
#include "../numerics/ParallelFor.hpp"
#include "../numerics/Random.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace PrecipitationStatistics {

enum class PitchDistribution { NarrowBin, UniformAngle, IsotropicHemisphere };

struct Interval {
    double low{};
    double high{};
};

inline Interval wilsonInterval(std::size_t successes, std::size_t trials) {
    if (trials == 0 || successes > trials)
        throw std::invalid_argument("Wilson interval requires 0 <= k <= N and N > 0");
    constexpr double z = 1.959963984540054;
    const double n = static_cast<double>(trials);
    const double fraction = static_cast<double>(successes) / n;
    const double denominator = 1.0 + z * z / n;
    const double center = (fraction + z * z / (2.0 * n)) / denominator;
    const double halfWidth = z / denominator *
        std::sqrt(fraction * (1.0 - fraction) / n + z * z / (4.0 * n * n));
    return {successes == 0 ? 0.0 : std::max(0.0, center - halfWidth),
            successes == trials ? 1.0 : std::min(1.0, center + halfWidth)};
}

struct Counts {
    std::size_t precipitated{};
    std::size_t mirrored{};
    std::size_t escaped{};
    std::size_t unresolved{};

    std::size_t total() const {
        return precipitated + mirrored + escaped + unresolved;
    }

    double precipitatingFraction() const {
        if (total() == 0) throw std::invalid_argument("Empty outcome sample");
        return static_cast<double>(precipitated) / static_cast<double>(total());
    }

    Interval precipitatingInterval() const {
        return wilsonInterval(precipitated, total());
    }
};

inline Counts count(const std::vector<SimulationResult>& results,
                    std::size_t limit) {
    if (limit == 0 || limit > results.size())
        throw std::invalid_argument("Invalid outcome prefix length");
    Counts counts;
    for (std::size_t index = 0; index < limit; ++index) {
        switch (results[index].outcome) {
            case ParticleOutcome::Precipitated: ++counts.precipitated; break;
            case ParticleOutcome::Mirrored: ++counts.mirrored; break;
            case ParticleOutcome::Escaped: ++counts.escaped; break;
            case ParticleOutcome::Unresolved: ++counts.unresolved; break;
        }
    }
    return counts;
}

struct Config {
    double colatitude_rad{};
    PitchDistribution distribution{PitchDistribution::NarrowBin};
    double binCenter_deg{};
    double binWidth_deg{};
    std::size_t samples{};
    std::uint64_t seed{};
    std::size_t threads{1};
    int stepsPerGyroperiod{50};
};

struct Result {
    double L{};
    double predictedCriticalPitch_deg{};
    double dt_s{};
    Counts counts;
    std::vector<InitialCondition> initialConditions;
    std::vector<SimulationResult> trajectories;
};

inline double theoreticalCriticalPitch_deg(double colatitude_rad,
                                            double atmosphereHeight_m = 150000.0) {
    if (!std::isfinite(colatitude_rad) || colatitude_rad <= 0.0 ||
        colatitude_rad >= std::numbers::pi / 2.0 ||
        !std::isfinite(atmosphereHeight_m) || atmosphereHeight_m < 0.0 ||
        atmosphereHeight_m >= 450000.0)
        throw std::invalid_argument("Invalid source or atmosphere altitude");
    const DipoleField field;
    const double sourceRadius = C::earthRadius_m + 450000.0;
    const double boundaryRadius = C::earthRadius_m + atmosphereHeight_m;
    const double sine = std::sin(colatitude_rad);
    const double L = sourceRadius / (C::earthRadius_m * sine * sine);
    const double boundaryColatitude = std::asin(std::sqrt(
        boundaryRadius / (L * C::earthRadius_m)));
    const Vector3 source(sourceRadius * sine, 0.0,
                         sourceRadius * std::cos(colatitude_rad));
    const Vector3 boundary(boundaryRadius * std::sin(boundaryColatitude), 0.0,
                           boundaryRadius * std::cos(boundaryColatitude));
    return std::asin(std::sqrt(field.evaluate(source).norm() /
                               field.evaluate(boundary).norm())) *
           180.0 / std::numbers::pi;
}

inline Result run(const Config& config) {
    if (!std::isfinite(config.colatitude_rad) || config.colatitude_rad <= 0.0 ||
        config.colatitude_rad >= std::numbers::pi / 2.0 ||
        config.samples == 0 || config.threads == 0 ||
        config.stepsPerGyroperiod < 4 ||
        (config.distribution == PitchDistribution::NarrowBin &&
         (!std::isfinite(config.binCenter_deg) ||
          !std::isfinite(config.binWidth_deg) || config.binWidth_deg <= 0.0 ||
          config.binCenter_deg - config.binWidth_deg / 2.0 <= 0.0 ||
          config.binCenter_deg + config.binWidth_deg / 2.0 >= 90.0))) {
        throw std::invalid_argument("PrecipitationStatistics: invalid configuration");
    }
    const DipoleField field;
    const double radius = C::earthRadius_m + 450000.0;
    const double sine = std::sin(config.colatitude_rad);
    const Vector3 position(radius * sine, 0.0,
                           radius * std::cos(config.colatitude_rad));
    const Vector3 magnetic = field.evaluate(position);
    const Vector3 along = magnetic * (-1.0 / magnetic.norm());
    const Vector3 across(0.0, 1.0, 0.0);
    const Vector3 around = along.cross(across);
    const double momentum = TrajectoryExperiments::momentumMagnitude();
    const double dt_s = TrajectoryExperiments::gyroPeriod(magnetic.norm()) /
        static_cast<double>(config.stepsPerGyroperiod);
    const SimulationEngine engine(field, -C::elementaryCharge_C, C::electronMass_kg,
        {dt_s, 150000.0, C::earthRadius_m + 800000.0, 0.05,
         static_cast<std::uint64_t>(std::ceil(0.05 / dt_s)) + 1});
    Result result;
    result.L = radius / (C::earthRadius_m * sine * sine);
    result.predictedCriticalPitch_deg = theoreticalCriticalPitch_deg(
        config.colatitude_rad);
    result.dt_s = dt_s;
    result.initialConditions.reserve(config.samples);
    Random random(config.seed);
    for (std::size_t index = 0; index < config.samples; ++index) {
        const double u = random.uniform01();
        double pitch_rad{};
        switch (config.distribution) {
            case PitchDistribution::NarrowBin:
                pitch_rad = (config.binCenter_deg +
                    (u - 0.5) * config.binWidth_deg) * std::numbers::pi / 180.0;
                break;
            case PitchDistribution::UniformAngle:
                pitch_rad = std::max(u, 1e-15) * std::numbers::pi / 2.0;
                break;
            case PitchDistribution::IsotropicHemisphere:
                // The same quantile u can be paired with UniformAngle:
                // cos(alpha) remains uniform, while both pitches grow with u.
                pitch_rad = std::acos(1.0 - std::max(u, 1e-15));
                break;
        }
        const double phase_rad = 2.0 * std::numbers::pi * random.uniform01();
        const Vector3 direction = along * std::cos(pitch_rad) +
            (across * std::cos(phase_rad) + around * std::sin(phase_rad)) *
            std::sin(pitch_rad);
        result.initialConditions.push_back({index + 1,
            {position, direction * momentum, 0.0},
            TrajectoryExperiments::energy_eV, pitch_rad, phase_rad, 1.0});
    }
    result.trajectories.resize(config.samples);
    parallelFor(config.samples, config.threads, [&](std::size_t index) {
        result.trajectories[index] = engine.run(result.initialConditions[index].state);
    });
    result.counts = count(result.trajectories, result.trajectories.size());
    return result;
}

} // namespace PrecipitationStatistics

#endif // AURORA_PRECIPITATIONSTATISTICS_HPP
