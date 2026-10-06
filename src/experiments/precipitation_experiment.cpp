#include "PrecipitationStatistics.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>
#include <thread>

namespace {

const char* outcomeName(ParticleOutcome outcome) {
    switch (outcome) {
        case ParticleOutcome::Precipitated: return "Precipitated";
        case ParticleOutcome::Mirrored: return "Mirrored";
        case ParticleOutcome::Escaped: return "Escaped";
        case ParticleOutcome::Unresolved: return "Unresolved";
    }
    throw std::logic_error("Unknown outcome");
}

void writeCounts(std::ofstream& output, const PrecipitationStatistics::Counts& counts) {
    const auto interval = counts.precipitatingInterval();
    output << counts.total() << ',' << counts.precipitated << ',' << counts.mirrored
           << ',' << counts.escaped << ',' << counts.unresolved << ','
           << counts.precipitatingFraction() << ',' << interval.low << ','
           << interval.high;
}

void writeParticles(std::ofstream& output, const std::string& group,
                    const PrecipitationStatistics::Result& run) {
    for (std::size_t index = 0; index < run.trajectories.size(); ++index) {
        const auto& initial = run.initialConditions[index];
        const auto& result = run.trajectories[index];
        output << group << ',' << initial.id << ','
               << initial.energy_eV << ','
               << initial.pitch_rad * 180.0 / std::numbers::pi << ','
               << initial.gyroPhase_rad * 180.0 / std::numbers::pi << ','
               << initial.state.time_s << ','
               << initial.state.r.x << ',' << initial.state.r.y << ','
               << initial.state.r.z << ','
               << initial.state.p.x << ',' << initial.state.p.y << ','
               << initial.state.p.z << ',' << initial.weight << ','
               << outcomeName(result.outcome) << ',' << result.eventState.time_s << ','
               << result.eventState.r.x << ',' << result.eventState.r.y << ','
               << result.eventState.r.z << ',' << result.steps << '\n';
    }
}

std::ofstream outputFile(const std::filesystem::path& directory,
                         const std::string& name) {
    std::ofstream output(directory / name);
    if (!output) throw std::runtime_error("Cannot write " + (directory / name).string());
    output << std::setprecision(17);
    return output;
}

} // namespace

int main(int argc, char* argv[]) {
    try {
        if (argc != 2 && argc != 4) {
            std::cerr << "Usage: aurora_precipitation_experiment OUTPUT_DIRECTORY "
                         "[--threads N]\n";
            return 2;
        }
        std::size_t threads = std::min<std::size_t>(
            8, std::max(1u, std::thread::hardware_concurrency()));
        if (argc == 4) {
            if (std::string(argv[2]) != "--threads")
                throw std::invalid_argument("Expected --threads N");
            threads = std::stoull(argv[3]);
            if (threads == 0 || threads > 64)
                throw std::invalid_argument("Thread count must be in [1,64]");
        }
        const std::filesystem::path directory(argv[1]);
        std::filesystem::create_directories(directory);
        auto bins = outputFile(directory, "precipitation_bins.csv");
        auto sources = outputFile(directory, "source_distributions.csv");
        auto convergence = outputFile(directory, "source_convergence.csv");
        auto particles = outputFile(directory, "precipitation_particles.csv");
        constexpr const char* countColumns =
            "N,precipitated,mirrored,escaped,unresolved,P_precip,"
            "wilson_low_95,wilson_high_95";
        bins << "group,colatitude_deg,L,pitch_center_deg,pitch_width_deg,"
                "theory_critical_deg,energy_eV,seed,dt_s," << countColumns << '\n';
        sources << "group,distribution,colatitude_deg,L,theory_critical_deg,"
                   "theory_fraction,energy_eV,seed,dt_s," << countColumns << '\n';
        convergence << "group,distribution,seed,prefix_N,precipitated,mirrored,"
                       "escaped,unresolved,P_precip,wilson_low_95,wilson_high_95\n";
        particles << "group,id,energy_eV,pitch_deg,phase_deg,initial_time_s,"
                     "initial_x_m,initial_y_m,initial_z_m,initial_px_kg_m_s,"
                     "initial_py_kg_m_s,initial_pz_kg_m_s,sample_weight,outcome,event_time_s,"
                     "event_x_m,event_y_m,event_z_m,steps\n";

        constexpr double colatitudes[]{0.25, 0.32, 0.40};
        constexpr double offsets_deg[]{-0.20, -0.10, 0.0, 0.10, 0.20};
        std::size_t groupIndex = 0;
        for (std::size_t line = 0; line < 3; ++line) {
            const double critical = PrecipitationStatistics::theoreticalCriticalPitch_deg(
                colatitudes[line]);
            for (std::size_t bin = 0; bin < 5; ++bin) {
                const std::string group = "bin_" + std::to_string(line) + "_" +
                                          std::to_string(bin);
                const std::uint64_t seed = 20261006 + groupIndex * 1009;
                const double center = critical + offsets_deg[bin];
                std::cout << "Computing " << group << ", center=" << center
                          << " deg" << std::endl;
                const PrecipitationStatistics::Config config{
                    colatitudes[line], PrecipitationStatistics::PitchDistribution::NarrowBin,
                    center, 0.10, 16, seed, threads, 50};
                const auto result = PrecipitationStatistics::run(config);
                bins << group << ',' << colatitudes[line] * 180.0 /
                        std::numbers::pi << ',' << result.L << ',' << center << ','
                     << config.binWidth_deg << ',' << result.predictedCriticalPitch_deg
                     << ',' << TrajectoryExperiments::energy_eV << ',' << seed << ','
                     << result.dt_s << ',';
                writeCounts(bins, result.counts);
                bins << '\n';
                writeParticles(particles, group, result);
                ++groupIndex;
            }
        }

        constexpr double sourceColatitude = 0.32;
        constexpr std::uint64_t seeds[]{20261006, 20261007};
        for (int distribution = 0; distribution < 2; ++distribution) {
            const auto kind = distribution == 0 ?
                PrecipitationStatistics::PitchDistribution::UniformAngle :
                PrecipitationStatistics::PitchDistribution::IsotropicHemisphere;
            const std::string name = distribution == 0 ?
                "uniform_angle" : "isotropic_hemisphere";
            for (int seedIndex = 0; seedIndex < 2; ++seedIndex) {
                const std::string group = name + "_" + std::to_string(seedIndex);
                const std::uint64_t seed = seeds[seedIndex];
                std::cout << "Computing " << group << ", N=128" << std::endl;
                const PrecipitationStatistics::Config config{
                    sourceColatitude, kind, 0.0, 0.0, 128, seed, threads, 50};
                const auto result = PrecipitationStatistics::run(config);
                const double critical_rad = result.predictedCriticalPitch_deg *
                                            std::numbers::pi / 180.0;
                const double theory = distribution == 0 ?
                    critical_rad / (std::numbers::pi / 2.0) :
                    1.0 - std::cos(critical_rad);
                sources << group << ',' << name << ',' << sourceColatitude *
                           180.0 / std::numbers::pi << ',' << result.L << ','
                        << result.predictedCriticalPitch_deg << ',' << theory << ','
                        << TrajectoryExperiments::energy_eV << ',' << seed << ','
                        << result.dt_s << ',';
                writeCounts(sources, result.counts);
                sources << '\n';
                writeParticles(particles, group, result);
                for (std::size_t prefix : {std::size_t{32}, std::size_t{64},
                                           std::size_t{128}}) {
                    const auto counts = PrecipitationStatistics::count(
                        result.trajectories, prefix);
                    const auto interval = counts.precipitatingInterval();
                    convergence << group << ',' << name << ',' << seed << ','
                                << prefix << ',' << counts.precipitated << ','
                                << counts.mirrored << ',' << counts.escaped << ','
                                << counts.unresolved << ','
                                << counts.precipitatingFraction() << ','
                                << interval.low << ',' << interval.high << '\n';
                }
            }
        }
        if (!bins || !sources || !convergence || !particles)
            throw std::runtime_error("Statistics CSV write failed");
        std::cout << "Saved statistical experiment CSV files to " << directory << '\n';
    } catch (const std::exception& error) {
        std::cerr << "Precipitation experiment failed: " << error.what() << '\n';
        return 1;
    }
}
