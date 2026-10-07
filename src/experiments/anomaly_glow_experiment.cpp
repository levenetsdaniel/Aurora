#include "PrecipitationStatistics.hpp"
#include "../visualization/AuroraPipeline.hpp"
#include "../io/AtmosphereCsvLoader.hpp"
#include "../io/IaaCrossSectionLoader.hpp"
#include "../io/KineticsCsvLoader.hpp"
#include "../io/SpectralLineCsvLoader.hpp"
#include "../physics/LocalizedDipoleField.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef AURORA_SOURCE_DIR
#error AURORA_SOURCE_DIR is required
#endif

namespace {

constexpr double colatitude_rad = 0.32;
constexpr double depth_m = 500000.0;
constexpr double surfaceFraction = 0.01;
constexpr double boundaryHeight_m = 150000.0;
constexpr double sourceHeight_m = 450000.0;
constexpr double sourceRate_s1 = 1.0;
constexpr std::array<std::uint64_t, 3> seeds{20261006, 20261007, 20261008};

const char* outcomeName(ParticleOutcome outcome) {
    switch (outcome) {
        case ParticleOutcome::Precipitated: return "Precipitated";
        case ParticleOutcome::Mirrored: return "Mirrored";
        case ParticleOutcome::Escaped: return "Escaped";
        case ParticleOutcome::Unresolved: return "Unresolved";
    }
    throw std::logic_error("Unknown outcome");
}

Vector3 sourcePoint(double altitude_m) {
    const double radius = C::earthRadius_m + altitude_m;
    return {radius * std::sin(colatitude_rad), 0.0,
            radius * std::cos(colatitude_rad)};
}

std::vector<SimulationResult> trace(const MagneticField& field,
                                    const std::vector<InitialCondition>& initial,
                                    double dt_s, std::size_t threads) {
    const SimulationEngine engine(field, -C::elementaryCharge_C, C::electronMass_kg,
        {dt_s, boundaryHeight_m, C::earthRadius_m + 800000.0, 0.05,
         static_cast<std::uint64_t>(std::ceil(0.05 / dt_s)) + 1});
    std::vector<SimulationResult> results(initial.size());
    parallelFor(initial.size(), threads, [&](std::size_t index) {
        results[index] = engine.run(initial[index].state);
    });
    return results;
}

std::vector<PrecipitatingElectron> precipitating(
    const std::vector<InitialCondition>& initial,
    const std::vector<SimulationResult>& results) {
    std::vector<PrecipitatingElectron> electrons;
    for (std::size_t index = 0; index < initial.size(); ++index) {
        if (results[index].outcome != ParticleOutcome::Precipitated) continue;
        const ParticleState& state = results[index].eventState;
        const Vector3 velocity = Relativity::velocityFromMomentum(
            state.p, C::electronMass_kg);
        const double cosine = -velocity.dot(state.r) /
            (velocity.norm() * state.r.norm());
        if (!std::isfinite(cosine) || cosine <= 0.0 || cosine > 1.0 + 1e-12)
            throw std::runtime_error("Invalid precipitation incidence angle");
        const double energy_eV = Relativity::kineticEnergyFromMomentum_J(
            state.p, C::electronMass_kg) / C::elementaryCharge_C;
        electrons.push_back({initial[index].id, state, energy_eV,
            std::clamp(cosine, 0.0, 1.0), sourceRate_s1 / initial.size()});
    }
    return electrons;
}

ExcitationGridConfig commonGrid(
    const std::vector<PrecipitatingElectron>& baseline,
    const std::vector<PrecipitatingElectron>& anomaly) {
    Vector3 lower = sourcePoint(boundaryHeight_m);
    Vector3 upper = lower;
    for (const auto* collection : {&baseline, &anomaly}) {
        for (const auto& electron : *collection) {
            const Vector3& point = electron.state.r;
            lower = {std::min(lower.x, point.x), std::min(lower.y, point.y),
                     std::min(lower.z, point.z)};
            upper = {std::max(upper.x, point.x), std::max(upper.y, point.y),
                     std::max(upper.z, point.z)};
        }
    }
    constexpr double horizontal_m = 50000.0;
    constexpr double vertical_m = 20000.0;
    return {{lower.x - 250000.0, lower.y - 250000.0, lower.z - 300000.0},
            {horizontal_m, horizontal_m, vertical_m},
            static_cast<std::size_t>(std::ceil((upper.x - lower.x + 500000.0) /
                                               horizontal_m)),
            static_cast<std::size_t>(std::ceil((upper.y - lower.y + 500000.0) /
                                               horizontal_m)),
            static_cast<std::size_t>(std::ceil((upper.z - lower.z + 450000.0) /
                                               vertical_m))};
}

void saveProfile(const VolumeEmissionGrid& emission, std::uint64_t seed,
                 const char* scenario, std::ofstream& profile,
                 std::ofstream& totals) {
    const ExcitationGridConfig& grid = emission.configuration();
    const double cellVolume_m3 = grid.cellSize_m.x * grid.cellSize_m.y *
                                 grid.cellSize_m.z;
    for (const auto& line : emission.lines()) {
        std::map<int, double> byAltitude;
        double total = 0.0;
        for (std::size_t z = 0; z < grid.cellsZ; ++z) {
            for (std::size_t y = 0; y < grid.cellsY; ++y) {
                for (std::size_t x = 0; x < grid.cellsX; ++x) {
                    const double rate = emission.emissivity(line.id, x, y, z) *
                                        cellVolume_m3;
                    if (rate == 0.0) continue;
                    const Vector3 center{
                        grid.origin_m.x + (x + 0.5) * grid.cellSize_m.x,
                        grid.origin_m.y + (y + 0.5) * grid.cellSize_m.y,
                        grid.origin_m.z + (z + 0.5) * grid.cellSize_m.z};
                    const double height_km = (center.norm() - C::earthRadius_m) / 1000.0;
                    const int bin_km = static_cast<int>(std::floor(height_km / 10.0)) * 10;
                    byAltitude[bin_km] += rate;
                    total += rate;
                }
            }
        }
        for (const auto& [height_km, rate] : byAltitude)
            profile << seed << ',' << scenario << ',' << line.id << ','
                    << height_km << ',' << rate << '\n';
        totals << seed << ',' << scenario << ',' << line.id << ',' << total << '\n';
    }
}

} // namespace

int main(int argc, char* argv[]) {
    try {
        if (argc != 3 && argc != 5) {
            std::cerr << "Usage: aurora_anomaly_glow_experiment OUTPUT_DIRECTORY "
                         "IAA_ROOT [--threads N]\n";
            return 2;
        }
        std::size_t threads = 4;
        if (argc == 5) {
            if (std::string(argv[3]) != "--threads")
                throw std::invalid_argument("Expected --threads N");
            threads = std::stoull(argv[4]);
            if (threads < 1 || threads > 64)
                throw std::invalid_argument("Thread count must be in [1,64]");
        }
        const std::filesystem::path output(argv[1]);
        std::filesystem::create_directories(output);
        std::ofstream particles(output / "anomaly_glow_particles.csv");
        std::ofstream summary(output / "anomaly_glow_summary.csv");
        std::ofstream profile(output / "anomaly_glow_profile.csv");
        std::ofstream totals(output / "anomaly_glow_totals.csv");
        std::ofstream diagnostics(output / "anomaly_glow_diagnostics.csv");
        if (!particles || !summary || !profile || !totals || !diagnostics)
            throw std::runtime_error("Cannot create anomaly-glow CSV files");
        for (auto* stream : {&particles, &summary, &profile, &totals, &diagnostics})
            *stream << std::setprecision(17);
        particles << "seed,id,pitch_deg,baseline_50,zero_50,weak_50,"
                     "baseline_100,weak_100,baseline_entry_x_m,baseline_entry_y_m,"
                     "baseline_entry_z_m,weak_entry_x_m,weak_entry_y_m,weak_entry_z_m\n";
        summary << "seed,scenario,steps_per_gyroperiod,N,precipitated,mirrored,"
                   "escaped,unresolved,precip_to_mirror,mirror_to_precip\n";
        profile << "seed,scenario,line,height_bin_km,photon_rate_s1\n";
        totals << "seed,scenario,line,photon_rate_s1\n";
        diagnostics << "seed,scenario,primaries,cascade_tracks,excitation_events,"
                       "initial_energy_eV_s1,energy_residual_eV_s1,"
                       "untracked_excitation_rate_s1,outside_grid_rate_s1\n";

        std::ofstream manifest(output / "anomaly_glow_manifest.txt");
        if (!manifest) throw std::runtime_error("Cannot create anomaly-glow manifest");
        manifest << "experiment=paired_anomaly_glow\n"
                 << "seeds=20261006,20261007,20261008\n"
                 << "samples_per_seed=64\n"
                 << "energy_eV=10000\n"
                 << "pitch_center_deg=69\n"
                 << "pitch_width_deg=1.2\n"
                 << "colatitude_rad=0.32\n"
                 << "source_altitude_m=450000\n"
                 << "handoff_altitude_m=150000\n"
                 << "anomaly_depth_m=500000\n"
                 << "anomaly_surface_fraction=0.01\n"
                 << "steps_per_gyroperiod=50,100\n"
                 << "emission_steps_per_gyroperiod=50\n"
                 << "source_rate_electrons_s1=1\n"
                 << "iaa_root=" << std::filesystem::absolute(argv[2]).string() << '\n';

        const std::filesystem::path root = AURORA_SOURCE_DIR;
        const AtmosphereModel atmosphere = AtmosphereCsvLoader::load(
            root / "data/atmosphere/nrlmsis20_2020-03-20_70N_0E_demo.csv");
        const CrossSectionDatabase crossSections = IaaCrossSectionLoader::load(argv[2], 25000.0);
        const KineticsCatalog states = KineticsCsvLoader::load(root / "data/kinetics");
        const SpectralLineCatalog lines = SpectralLineCsvLoader::load(
            root / "data/kinetics/lines.csv", states);
        const DipoleField background;
        const double surfacePeak_T = surfaceFraction *
            background.evaluate(sourcePoint(0.0)).norm();
        const LocalizedDipoleField zero(colatitude_rad, depth_m, 0.0);
        const LocalizedDipoleField weak(colatitude_rad, depth_m, surfacePeak_T);
        constexpr std::size_t samples = 64;
        constexpr double center_deg = 69.0;
        constexpr double width_deg = 1.2;

        for (const std::uint64_t seed : seeds) {
            std::cout << "Seed " << seed << ": tracing paired particles" << std::endl;
            const PrecipitationStatistics::Config config{
                colatitude_rad, PrecipitationStatistics::PitchDistribution::NarrowBin,
                center_deg, width_deg, samples, seed, threads, 50};
            const auto base = PrecipitationStatistics::run(config);
            const auto& initial = base.initialConditions;
            const auto baseFine = trace(background, initial, base.dt_s / 2.0, threads);
            const auto zeroCoarse = trace(zero, initial, base.dt_s, threads);
            const auto weakCoarse = trace(weak, initial, base.dt_s, threads);
            const auto weakFine = trace(weak, initial, base.dt_s / 2.0, threads);
            for (std::size_t index = 0; index < samples; ++index) {
                const auto& first = base.trajectories[index];
                const auto& control = zeroCoarse[index];
                if (first.outcome != control.outcome ||
                    (first.eventState.r - control.eventState.r).norm() != 0.0 ||
                    (first.eventState.p - control.eventState.p).norm() != 0.0)
                    throw std::runtime_error("Zero-anomaly paired control differs");
                const auto& changed = weakCoarse[index];
                particles << seed << ',' << initial[index].id << ','
                          << initial[index].pitch_rad * 180.0 / std::numbers::pi << ','
                          << outcomeName(first.outcome) << ',' << outcomeName(control.outcome)
                          << ',' << outcomeName(changed.outcome) << ','
                          << outcomeName(baseFine[index].outcome) << ','
                          << outcomeName(weakFine[index].outcome) << ','
                          << first.eventState.r.x << ',' << first.eventState.r.y << ','
                          << first.eventState.r.z << ',' << changed.eventState.r.x << ','
                          << changed.eventState.r.y << ',' << changed.eventState.r.z << '\n';
            }
            const auto saveCounts = [&](const char* name, int steps,
                                        const std::vector<SimulationResult>& results,
                                        const std::vector<SimulationResult>& reference) {
                const auto counts = PrecipitationStatistics::count(results, samples);
                std::size_t precipToMirror = 0;
                std::size_t mirrorToPrecip = 0;
                for (std::size_t index = 0; index < samples; ++index) {
                    if (reference[index].outcome == ParticleOutcome::Precipitated &&
                        results[index].outcome == ParticleOutcome::Mirrored) ++precipToMirror;
                    if (reference[index].outcome == ParticleOutcome::Mirrored &&
                        results[index].outcome == ParticleOutcome::Precipitated) ++mirrorToPrecip;
                }
                summary << seed << ',' << name << ',' << steps << ',' << samples << ','
                        << counts.precipitated << ',' << counts.mirrored << ','
                        << counts.escaped << ',' << counts.unresolved << ','
                        << precipToMirror << ',' << mirrorToPrecip << '\n';
            };
            saveCounts("baseline", 50, base.trajectories, base.trajectories);
            saveCounts("zero", 50, zeroCoarse, base.trajectories);
            saveCounts("weak", 50, weakCoarse, base.trajectories);
            saveCounts("baseline", 100, baseFine, baseFine);
            saveCounts("weak", 100, weakFine, baseFine);

            const auto baseElectrons = precipitating(initial, base.trajectories);
            const auto weakElectrons = precipitating(initial, weakCoarse);
            const ExcitationGridConfig grid = commonGrid(baseElectrons, weakElectrons);
            const AuroraPipeline pipeline(atmosphere, crossSections, states, lines,
                {grid, {120000.0, 250.0, 20000, 2048, 1.0}, threads, 64, seed});
            for (const auto& [name, electrons] :
                 {std::pair<const char*, const std::vector<PrecipitatingElectron>&>{
                      "baseline", baseElectrons},
                  {"weak", weakElectrons}}) {
                std::cout << "Seed " << seed << ": atmosphere " << name << " ("
                          << electrons.size() << " primaries)" << std::endl;
                if (electrons.empty())
                    throw std::runtime_error("No precipitating electrons for emission");
                const auto result = pipeline.run(electrons);
                if (result.outsideGridRate_s1 > 0.0)
                    throw std::runtime_error("Emission grid lost excitation events");
                diagnostics << seed << ',' << name << ',' << result.primaryCount << ','
                            << result.cascadeTracks << ',' << result.excitationEvents << ','
                            << result.initialEnergy_eV_s1 << ','
                            << result.energyResidual_eV_s1 << ','
                            << result.untrackedExcitationRate_s1 << ','
                            << result.outsideGridRate_s1 << '\n';
                result.emission.saveSparseCsv(output /
                    ("anomaly_glow_" + std::to_string(seed) + "_" + name + ".csv"));
                saveProfile(result.emission, seed, name, profile, totals);
            }
            for (auto* stream : {&particles, &summary, &profile, &totals, &diagnostics})
                stream->flush();
        }
        if (!particles || !summary || !profile || !totals || !diagnostics || !manifest)
            throw std::runtime_error("Anomaly-glow CSV write failed");
        std::cout << "Saved anomaly-glow results to " << output << '\n';
    } catch (const std::exception& error) {
        std::cerr << "Anomaly-glow experiment failed: " << error.what() << '\n';
        return 1;
    }
}
