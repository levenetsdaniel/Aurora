#include "ElectronSceneData.hpp"
#include "VtkElectronScene.hpp"
#include "AuroraPipeline.hpp"
#include "EmissionCsvReplay.hpp"
#include "ProgressReporter.hpp"
#include "../ensemble/ArcElectronSource.hpp"
#include "../io/AtmosphereCsvLoader.hpp"
#include "../io/IaaCrossSectionLoader.hpp"
#include "../io/KineticsCsvLoader.hpp"
#include "../io/SpectralLineCsvLoader.hpp"
#include "../emission/OpticalRenderer.hpp"
#include "../numerics/Random.hpp"
#include "../numerics/ParallelFor.hpp"
#include "../physics/DipoleField.hpp"
#include "../physics/LocalizedDipoleField.hpp"
#include "../physics/Relativity.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>

#ifndef AURORA_SOURCE_DIR
#error AURORA_SOURCE_DIR is required
#endif

namespace {

const char* outcomeName(ParticleOutcome outcome) {
    switch (outcome) {
        case ParticleOutcome::Precipitated: return "Precipitated";
        case ParticleOutcome::Mirrored: return "Mirrored";
        case ParticleOutcome::Escaped: return "Escaped";
        case ParticleOutcome::Unresolved: return "Unresolved";
    }
    throw std::logic_error("Unknown particle outcome");
}

void saveParticleOutcomes(const ElectronSceneData& scene,
                          const std::filesystem::path& path) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("Cannot write particle outcomes");
    output << "id,energy_eV,pitch_deg,phase_deg,initial_time_s,initial_x_m,initial_y_m,initial_z_m,"
              "initial_px_kg_m_s,initial_py_kg_m_s,initial_pz_kg_m_s,outcome,"
              "event_time_s,event_x_m,event_y_m,event_z_m,weight_electrons_s1\n";
    output << std::setprecision(17);
    for (const ElectronPath& particle : scene.paths) {
        output << particle.id << ',' << particle.energy_eV << ','
               << particle.initialPitch_rad * 180.0 / std::numbers::pi << ','
               << particle.initialGyroPhase_rad * 180.0 / std::numbers::pi << ','
               << particle.initialState.time_s << ','
               << particle.initialState.r.x << ','
               << particle.initialState.r.y << ','
               << particle.initialState.r.z << ','
               << particle.initialState.p.x << ','
               << particle.initialState.p.y << ','
               << particle.initialState.p.z << ','
               << outcomeName(particle.outcome) << ','
               << particle.terminalState.time_s << ','
               << particle.terminalState.r.x << ','
               << particle.terminalState.r.y << ','
               << particle.terminalState.r.z << ',' << particle.weight_s1 << '\n';
    }
    if (!output) throw std::runtime_error("Particle outcome CSV write failed");
}

std::vector<InitialCondition> makeInitialConditions(
    const MagneticField& sourceField, std::size_t markerCount,
    std::uint64_t seed, double sourceRate_s1, bool auroralArc,
    double minPitch_rad, double maxPitch_rad) {
    Random random(seed);
    const ArcElectronSource source(sourceField, auroralArc, minPitch_rad,
                                   maxPitch_rad, sourceRate_s1, markerCount);
    std::vector<InitialCondition> initialConditions;
    initialConditions.reserve(markerCount);
    for (std::size_t index = 0; index < markerCount; ++index)
        initialConditions.push_back(source.sample(random, index + 1));
    return initialConditions;
}

ElectronSceneData traceScene(const MagneticField& field,
                             const std::vector<InitialCondition>& initialConditions,
                             std::size_t threads, const std::string& stage) {
    constexpr double atmosphereHeight_m = 150000.0;
    const SimulationConfig config{2e-8, atmosphereHeight_m,
        C::earthRadius_m + 800000.0, 0.05, 2500000};
    const SimulationEngine engine(field, -C::elementaryCharge_C,
                                  C::electronMass_kg, config);
    const ElectronSceneSampler sampler(engine, field);
    ElectronSceneData scene;
    scene.atmosphereHeight_m = atmosphereHeight_m;
    scene.paths.resize(initialConditions.size());
    ProgressReporter orbitProgress(stage, initialConditions.size(), std::cerr);
    parallelFor(initialConditions.size(), threads, [&](std::size_t index) {
        scene.paths[index] = sampler.traceElectron(initialConditions[index], 250);
        orbitProgress.advance();
    });
    for (std::size_t index = 0; index < initialConditions.size(); ++index) {
        if (index % std::max<std::size_t>(1, initialConditions.size() / 12) == 0) {
            scene.magneticLines_m.push_back(
                sampler.traceMagneticLine(initialConditions[index].state.r));
        }
    }
    return scene;
}

std::vector<PrecipitatingElectron> precipitation(const ElectronSceneData& scene) {
    std::vector<PrecipitatingElectron> electrons;
    for (const ElectronPath& path : scene.paths) {
        if (path.outcome != ParticleOutcome::Precipitated) continue;
        const ParticleState& state = path.terminalState;
        const Vector3 velocity = Relativity::velocityFromMomentum(
            state.p, C::electronMass_kg);
        const double cosine = -velocity.dot(state.r) /
            (velocity.norm() * state.r.norm());
        const double energy_eV = Relativity::kineticEnergyFromMomentum_J(
            state.p, C::electronMass_kg) / C::elementaryCharge_C;
        if (!std::isfinite(cosine) || cosine <= 0.0 || cosine > 1.0 + 1e-12) {
            throw std::runtime_error("Precipitating electron has invalid incidence angle");
        }
        electrons.push_back({path.id, state, energy_eV,
                             std::clamp(cosine, 0.0, 1.0), path.weight_s1});
    }
    return electrons;
}

ExcitationGridConfig emissionGrid(const std::vector<PrecipitatingElectron>& electrons) {
    if (electrons.empty()) throw std::invalid_argument("No precipitating electrons");
    Vector3 lower = electrons.front().state.r;
    Vector3 upper = lower;
    for (const auto& electron : electrons) {
        const Vector3& point = electron.state.r;
        lower = {std::min(lower.x, point.x), std::min(lower.y, point.y),
                 std::min(lower.z, point.z)};
        upper = {std::max(upper.x, point.x), std::max(upper.y, point.y),
                 std::max(upper.z, point.z)};
    }
    constexpr double horizontalCell_m = 60000.0;
    constexpr double verticalCell_m = 30000.0;
    const Vector3 origin{lower.x - 250000.0, lower.y - 250000.0,
                         lower.z - 200000.0};
    return {origin, {horizontalCell_m, horizontalCell_m, verticalCell_m},
        static_cast<std::size_t>(std::ceil((upper.x - lower.x + 500000.0) /
                                           horizontalCell_m)),
        static_cast<std::size_t>(std::ceil((upper.y - lower.y + 500000.0) /
                                           horizontalCell_m)),
        static_cast<std::size_t>(std::ceil((upper.z - lower.z + 300000.0) /
                                           verticalCell_m))};
}

AuroraPipelineResult makeAurora(const ElectronSceneData& scene,
                                const std::filesystem::path& iaaRoot,
                                std::size_t threads, std::uint64_t seed,
                                std::size_t cascadeMarkers,
                                const std::vector<std::uint64_t>& pairedSourceIds = {}) {
    const std::filesystem::path sourceRoot = AURORA_SOURCE_DIR;
    const AtmosphereModel atmosphere = AtmosphereCsvLoader::load(
        sourceRoot / "data/atmosphere/nrlmsis20_2020-03-20_70N_0E_demo.csv");
    const CrossSectionDatabase crossSections = IaaCrossSectionLoader::load(
        iaaRoot, 25000.0);
    const KineticsCatalog states = KineticsCsvLoader::load(
        sourceRoot / "data/kinetics");
    const SpectralLineCatalog lines = SpectralLineCsvLoader::load(
        sourceRoot / "data/kinetics/lines.csv", states);
    auto electrons = precipitation(scene);
    if (!pairedSourceIds.empty()) {
        const std::unordered_set<std::uint64_t> chosen(
            pairedSourceIds.begin(), pairedSourceIds.end());
        std::vector<PrecipitatingElectron> selected;
        selected.reserve(pairedSourceIds.size());
        const double scale = static_cast<double>(scene.paths.size()) /
                             static_cast<double>(pairedSourceIds.size());
        for (auto& electron : electrons) {
            if (!chosen.contains(electron.id)) continue;
            electron.weight *= scale;
            selected.push_back(electron);
        }
        electrons = std::move(selected);
    } else if (cascadeMarkers > 0 && cascadeMarkers < electrons.size()) {
        std::vector<PrecipitatingElectron> selected;
        selected.reserve(cascadeMarkers);
        double totalRate_s1 = 0.0;
        for (const auto& electron : electrons) totalRate_s1 += electron.weight;
        for (std::size_t index = 0; index < cascadeMarkers; ++index) {
            const std::size_t source = (2 * index + 1) * electrons.size() /
                                       (2 * cascadeMarkers);
            selected.push_back(electrons[source]);
            selected.back().weight = totalRate_s1 /
                                     static_cast<double>(cascadeMarkers);
        }
        electrons = std::move(selected);
    }
    const AuroraPipelineConfig config{
        emissionGrid(electrons),
        {120000.0, 250.0, 20000, 2048, 1.0},
        threads, 64, seed};
    const AuroraPipeline pipeline(atmosphere, crossSections, states, lines, config);
    ProgressReporter cascadeProgress("Atmospheric cascades", electrons.size(), std::cerr);
    return pipeline.run(electrons, [&] { cascadeProgress.advance(); });
}

std::vector<std::uint64_t> pairedSourceIds(std::size_t markerCount,
                                            std::size_t cascadeMarkers) {
    if (cascadeMarkers == 0 || cascadeMarkers >= markerCount) return {};
    std::vector<std::uint64_t> ids;
    ids.reserve(cascadeMarkers);
    for (std::size_t index = 0; index < cascadeMarkers; ++index)
        ids.push_back(1 + (2 * index + 1) * markerCount / (2 * cascadeMarkers));
    return ids;
}

void saveSpectralViews(const ElectronSceneData& scene,
                       const VolumeEmissionGrid& emission,
                       const std::filesystem::path& outputDirectory) {
    Vector3 center;
    std::size_t impacts = 0;
    for (const ElectronPath& path : scene.paths) {
        if (path.outcome != ParticleOutcome::Precipitated) continue;
        center = center + path.terminalState.r;
        ++impacts;
    }
    if (impacts == 0) return;
    center = center / static_cast<double>(impacts);
    const Vector3 radial = center / center.norm();
    const OpticalRenderer camera({center + radial * 1500000.0,
                                  radial * -1.0, {0.0, 1.0, 0.0},
                                  60.0, 320, 200});
    const auto images = camera.render(emission);
    OpticalRenderer::saveSpectralCsv(images, outputDirectory);
    const auto findLine = [&](const std::string& id) -> const SpectralImage& {
        const auto found = std::find_if(images.begin(), images.end(),
            [&](const SpectralImage& image) { return image.lineId == id; });
        if (found == images.end()) throw std::runtime_error("Missing spectral image " + id);
        return *found;
    };
    const SpectralImage& red = findLine("O_6300");
    const SpectralImage& green = findLine("O_5577");
    double maximum = 0.0;
    for (double value : red.radiance_ph_m2_s1_sr1) maximum = std::max(maximum, value);
    for (double value : green.radiance_ph_m2_s1_sr1) maximum = std::max(maximum, value);
    const auto preview = OpticalRenderer::makeRgbPreview(
        red, green, maximum > 0.0 ? maximum / 2.0 : 1.0);
    preview.savePpm(outputDirectory / "aurora_preview.ppm");
}

double emissionPeak(const VolumeEmissionGrid& emission) {
    double greenPeak = 0.0;
    double redPeak = 0.0;
    const auto& grid = emission.configuration();
    for (std::size_t z = 0; z < grid.cellsZ; ++z) {
        for (std::size_t y = 0; y < grid.cellsY; ++y) {
            for (std::size_t x = 0; x < grid.cellsX; ++x) {
                greenPeak = std::max(greenPeak,
                    emission.emissivity("O_5577", x, y, z));
                redPeak = std::max(redPeak,
                    emission.emissivity("O_6300", x, y, z));
            }
        }
    }
    return std::max(greenPeak, redPeak);
}

void addGlow(VtkElectronScene& visualizer, const VolumeEmissionGrid& emission,
             double sharedReferencePeak = 0.0) {
    double greenPeak = 0.0;
    double redPeak = 0.0;
    const auto& grid = emission.configuration();
    for (std::size_t z = 0; z < grid.cellsZ; ++z) {
        for (std::size_t y = 0; y < grid.cellsY; ++y) {
            for (std::size_t x = 0; x < grid.cellsX; ++x) {
                greenPeak = std::max(greenPeak,
                    emission.emissivity("O_5577", x, y, z));
                redPeak = std::max(redPeak,
                    emission.emissivity("O_6300", x, y, z));
            }
        }
    }
    const double referencePeak = sharedReferencePeak > 0.0 ?
        sharedReferencePeak : std::max(greenPeak, redPeak);
    const double greenOpacity = referencePeak > 0.0 ?
        std::sqrt(greenPeak / referencePeak) : 0.0;
    const double redOpacity = referencePeak > 0.0 ?
        std::sqrt(redPeak / referencePeak) : 0.0;
    visualizer.setCameraScale(0.31);
    visualizer.addEmission(emission, "O_6300", 1.0, 0.16, 0.10, redOpacity);
    visualizer.addEmission(emission, "O_5577", 0.14, 1.0, 0.27, greenOpacity);
}

bool hasSelectedPrecipitation(const ElectronSceneData& scene,
                              const std::vector<std::uint64_t>& ids) {
    if (ids.empty()) return scene.precipitationCount() > 0;
    const std::unordered_set<std::uint64_t> chosen(ids.begin(), ids.end());
    return std::any_of(scene.paths.begin(), scene.paths.end(),
        [&](const ElectronPath& path) {
            return path.outcome == ParticleOutcome::Precipitated &&
                   chosen.contains(path.id);
        });
}

void savePairedOutcomes(const ElectronSceneData& baseline,
                        const ElectronSceneData& anomaly,
                        const std::filesystem::path& path) {
    if (baseline.paths.size() != anomaly.paths.size())
        throw std::invalid_argument("Paired scenes have different marker counts");
    std::ofstream output(path);
    if (!output) throw std::runtime_error("Cannot write paired outcomes");
    output << std::setprecision(17);
    output << "id,initial_x_m,initial_y_m,initial_z_m,initial_px_kg_m_s,"
              "initial_py_kg_m_s,initial_pz_kg_m_s,weight_electrons_s1,"
              "baseline_outcome,anomaly_outcome,baseline_event_time_s,"
              "anomaly_event_time_s,baseline_event_x_m,baseline_event_y_m,"
              "baseline_event_z_m,anomaly_event_x_m,anomaly_event_y_m,"
              "anomaly_event_z_m,paired_entry_shift_m\n";
    for (std::size_t index = 0; index < baseline.paths.size(); ++index) {
        const ElectronPath& first = baseline.paths[index];
        const ElectronPath& second = anomaly.paths[index];
        if (first.id != second.id || first.weight_s1 != second.weight_s1 ||
            first.initialState.time_s != second.initialState.time_s ||
            (first.initialState.r - second.initialState.r).norm() != 0.0 ||
            (first.initialState.p - second.initialState.p).norm() != 0.0)
            throw std::runtime_error("Paired scenes changed an initial state");
        const bool bothEntered = first.outcome == ParticleOutcome::Precipitated &&
                                 second.outcome == ParticleOutcome::Precipitated;
        output << first.id << ',' << first.initialState.r.x << ','
               << first.initialState.r.y << ',' << first.initialState.r.z << ','
               << first.initialState.p.x << ',' << first.initialState.p.y << ','
               << first.initialState.p.z << ',' << first.weight_s1 << ','
               << outcomeName(first.outcome) << ',' << outcomeName(second.outcome)
               << ',' << first.terminalState.time_s << ','
               << second.terminalState.time_s << ','
               << first.terminalState.r.x << ',' << first.terminalState.r.y << ','
               << first.terminalState.r.z << ',' << second.terminalState.r.x
               << ',' << second.terminalState.r.y << ',' << second.terminalState.r.z
               << ',';
        if (bothEntered)
            output << (first.terminalState.r - second.terminalState.r).norm();
        output << '\n';
    }
    if (!output) throw std::runtime_error("Paired outcome CSV write failed");
}

void saveFrames(const ElectronSceneData& scene, VtkElectronScene& visualizer,
                const std::filesystem::path& directory, std::size_t frameCount,
                double lastTime_s, const std::string& stage) {
    std::filesystem::create_directories(directory);
    ProgressReporter frameProgress(stage, frameCount, std::cerr);
    for (std::size_t frame = 0; frame < frameCount; ++frame) {
        const double time_s = ElectronAnimation::frameTime_s(
            frame, frameCount, lastTime_s);
        visualizer.setAnimationFrame(scene, time_s);
        std::ostringstream name;
        name << "frame_" << std::setw(4) << std::setfill('0') << frame << ".png";
        visualizer.savePng(directory / name.str());
        frameProgress.advance();
    }
}

} // namespace

int main(int argc, char** argv) {
    try {
        std::filesystem::path output = "electron_flow.png";
        std::filesystem::path framesDirectory;
        std::size_t frameCount = 80;
        int width = 1600;
        int height = 1000;
        std::size_t markerCount = 48;
        std::size_t cascadeMarkers = 0;
        std::uint64_t seed = 20260928;
        std::size_t threads = std::max(1u, std::thread::hardware_concurrency());
        double sourceRate_s1 = 1.0e20;
        double minPitch_deg = 0.0;
        double maxPitch_deg = 90.0;
        std::filesystem::path iaaRoot;
        std::filesystem::path emissionCsv;
        bool aurora = false;
        bool compareAnomaly = false;
        double anomalyDepth_km = 500.0;
        double anomalySurfaceFraction = 0.01;
        bool interactive = false;
        for (int index = 1; index < argc; ++index) {
            const std::string argument = argv[index];
            if (argument == "--output" && index + 1 < argc) {
                output = argv[++index];
            } else if (argument == "--frames-dir" && index + 1 < argc) {
                framesDirectory = argv[++index];
            } else if (argument == "--frames" && index + 1 < argc) {
                frameCount = std::stoull(argv[++index]);
            } else if (argument == "--width" && index + 1 < argc) {
                width = std::stoi(argv[++index]);
            } else if (argument == "--height" && index + 1 < argc) {
                height = std::stoi(argv[++index]);
            } else if (argument == "--markers" && index + 1 < argc) {
                markerCount = std::stoull(argv[++index]);
            } else if (argument == "--seed" && index + 1 < argc) {
                seed = std::stoull(argv[++index]);
            } else if (argument == "--threads" && index + 1 < argc) {
                threads = std::stoull(argv[++index]);
            } else if (argument == "--cascade-markers" && index + 1 < argc) {
                cascadeMarkers = std::stoull(argv[++index]);
            } else if (argument == "--source-rate" && index + 1 < argc) {
                sourceRate_s1 = std::stod(argv[++index]);
            } else if (argument == "--pitch-min-deg" && index + 1 < argc) {
                minPitch_deg = std::stod(argv[++index]);
            } else if (argument == "--pitch-max-deg" && index + 1 < argc) {
                maxPitch_deg = std::stod(argv[++index]);
            } else if (argument == "--iaa-root" && index + 1 < argc) {
                iaaRoot = argv[++index];
            } else if (argument == "--emission-csv" && index + 1 < argc) {
                emissionCsv = argv[++index];
            } else if (argument == "--compare-anomaly") {
                compareAnomaly = true;
            } else if (argument == "--anomaly-depth-km" && index + 1 < argc) {
                anomalyDepth_km = std::stod(argv[++index]);
            } else if (argument == "--anomaly-surface-fraction" && index + 1 < argc) {
                anomalySurfaceFraction = std::stod(argv[++index]);
            } else if (argument == "--aurora") {
                aurora = true;
            } else if (argument == "--interactive") {
                interactive = true;
            } else {
                throw std::invalid_argument(
                    "Usage: aurora_vtk_visualizer [--output image.png] "
                    "[--markers N] [--cascade-markers N] [--seed N] "
                    "[--threads N] [--source-rate N] "
                    "[--pitch-min-deg N] [--pitch-max-deg N] "
                    "[--aurora --iaa-root DIR | --emission-csv FILE] [--interactive] "
                    "[--compare-anomaly --anomaly-depth-km KM "
                    "--anomaly-surface-fraction F] "
                    "[--frames-dir EMPTY_DIR --frames N --width PX --height PX]");
            }
        }
        if (markerCount == 0 || markerCount > 5000 ||
            threads == 0 || threads > 64 ||
            !std::isfinite(sourceRate_s1) || sourceRate_s1 <= 0.0 ||
            !std::isfinite(minPitch_deg) || !std::isfinite(maxPitch_deg) ||
            minPitch_deg < 0.0 || minPitch_deg >= maxPitch_deg ||
            maxPitch_deg > 90.0 ||
            frameCount < 2 || frameCount > 600 ||
            width < 320 || width > 3840 || height < 240 || height > 2160 ||
            output.empty() || output.extension() != ".png" ||
            (aurora && iaaRoot.empty()) ||
            (aurora && !emissionCsv.empty()) ||
            (compareAnomaly && (!emissionCsv.empty() || interactive)) ||
            !std::isfinite(anomalyDepth_km) || anomalyDepth_km <= 0.0 ||
            anomalyDepth_km * 1000.0 >= C::earthRadius_m ||
            !std::isfinite(anomalySurfaceFraction) ||
            anomalySurfaceFraction < 0.0 || anomalySurfaceFraction > 1.0) {
            throw std::invalid_argument("Invalid marker count, rate, thread count, IAA root or output path");
        }
        if (!framesDirectory.empty() &&
            std::filesystem::exists(framesDirectory) &&
            !std::filesystem::is_empty(framesDirectory))
            throw std::invalid_argument("--frames-dir must be empty");
        if (!emissionCsv.empty())
            EmissionCsvReplay::validateCompanionManifest(
                emissionCsv, seed, markerCount, sourceRate_s1,
                minPitch_deg, maxPitch_deg);
        const DipoleField baseField;
        const auto initialConditions = makeInitialConditions(
            baseField, markerCount, seed, sourceRate_s1,
            aurora || compareAnomaly || !emissionCsv.empty(),
            minPitch_deg * std::numbers::pi / 180.0,
            maxPitch_deg * std::numbers::pi / 180.0);
        if (compareAnomaly) {
            // The entire source is generated once in the background field.
            // Both orbit integrators receive identical phase-space states.
            const double colatitude = 0.32;
            const Vector3 surface(
                C::earthRadius_m * std::sin(colatitude), 0.0,
                C::earthRadius_m * std::cos(colatitude));
            const double peak_T = anomalySurfaceFraction *
                                  baseField.evaluate(surface).norm();
            const LocalizedDipoleField anomalyField(
                colatitude, anomalyDepth_km * 1000.0, peak_T);
            const auto baseline = traceScene(
                baseField, initialConditions, threads, "Base trajectories");
            const auto changed = traceScene(
                anomalyField, initialConditions, threads, "Anomaly trajectories");
            if (!output.parent_path().empty())
                std::filesystem::create_directories(output.parent_path());
            const std::string stem = output.stem().string();
            const auto basePath = output.parent_path() / (stem + "_base.png");
            const auto changedPath = output.parent_path() / (stem + "_anomaly.png");
            const auto comparisonPath = output.parent_path() /
                                        (stem + "_comparison.csv");
            saveParticleOutcomes(baseline, output.parent_path() /
                (stem + "_base_particles.csv"));
            saveParticleOutcomes(changed, output.parent_path() /
                (stem + "_anomaly_particles.csv"));
            savePairedOutcomes(baseline, changed, comparisonPath);

            const auto selectedIds = pairedSourceIds(markerCount, cascadeMarkers);
            std::optional<AuroraPipelineResult> baseEmission;
            std::optional<AuroraPipelineResult> changedEmission;
            if (aurora && hasSelectedPrecipitation(baseline, selectedIds))
                baseEmission.emplace(makeAurora(
                    baseline, iaaRoot, threads, seed, 0, selectedIds));
            if (aurora && hasSelectedPrecipitation(changed, selectedIds))
                changedEmission.emplace(makeAurora(
                    changed, iaaRoot, threads, seed, 0, selectedIds));
            const double sharedPeak = std::max(
                baseEmission ? emissionPeak(baseEmission->emission) : 0.0,
                changedEmission ? emissionPeak(changedEmission->emission) : 0.0);
            VtkElectronScene baseView(baseline, width, height,
                                      !framesDirectory.empty());
            VtkElectronScene changedView(changed, width, height,
                                         !framesDirectory.empty());
            if (aurora) {
                baseView.setCameraScale(0.31);
                changedView.setCameraScale(0.31);
            }
            if (baseEmission) {
                addGlow(baseView, baseEmission->emission, sharedPeak);
                baseEmission->emission.saveSparseCsv(output.parent_path() /
                    (stem + "_base_emission.csv"));
                saveSpectralViews(baseline, baseEmission->emission,
                    output.parent_path() / (stem + "_base_spectral"));
            }
            if (changedEmission) {
                addGlow(changedView, changedEmission->emission, sharedPeak);
                changedEmission->emission.saveSparseCsv(output.parent_path() /
                    (stem + "_anomaly_emission.csv"));
                saveSpectralViews(changed, changedEmission->emission,
                    output.parent_path() / (stem + "_anomaly_spectral"));
            }
            if (!framesDirectory.empty()) {
                const double lastTime_s = std::max(
                    ElectronAnimation::lastEventTime_s(baseline),
                    ElectronAnimation::lastEventTime_s(changed));
                saveFrames(baseline, baseView, framesDirectory / "base",
                           frameCount, lastTime_s, "Base animation frames");
                saveFrames(changed, changedView, framesDirectory / "anomaly",
                           frameCount, lastTime_s, "Anomaly animation frames");
            }
            baseView.savePng(basePath);
            changedView.savePng(changedPath);
            const auto manifestPath = output.parent_path() /
                                      (stem + "_comparison_manifest.txt");
            std::ofstream manifest(manifestPath);
            if (!manifest) throw std::runtime_error("Cannot write comparison manifest");
            manifest << std::setprecision(17)
                     << "seed=" << seed << '\n'
                     << "visual_markers=" << markerCount << '\n'
                     << "source_rate_electrons_s1=" << sourceRate_s1 << '\n'
                     << "pitch_min_deg=" << minPitch_deg << '\n'
                     << "pitch_max_deg=" << maxPitch_deg << '\n'
                     << "anomaly_depth_km=" << anomalyDepth_km << '\n'
                     << "anomaly_surface_fraction=" << anomalySurfaceFraction << '\n'
                     << "anomaly_colatitude_rad=" << colatitude << '\n'
                     << "atmosphere_height_m=150000\n"
                     << "orbit_dt_s=2e-8\n"
                     << "shared_source_ids=" << (selectedIds.empty() ? markerCount :
                                                 selectedIds.size()) << '\n'
                     << "baseline_precipitated=" << baseline.precipitationCount() << '\n'
                     << "anomaly_precipitated=" << changed.precipitationCount() << '\n'
                     << "baseline_mirrored="
                     << baseline.countOutcome(ParticleOutcome::Mirrored) << '\n'
                     << "anomaly_mirrored="
                     << changed.countOutcome(ParticleOutcome::Mirrored) << '\n'
                     << "baseline_emission_recomputed=" << (baseEmission ? 1 : 0) << '\n'
                     << "anomaly_emission_recomputed=" << (changedEmission ? 1 : 0) << '\n'
                     << "baseline_cascade_primaries="
                     << (baseEmission ? baseEmission->primaryCount : 0) << '\n'
                     << "anomaly_cascade_primaries="
                     << (changedEmission ? changedEmission->primaryCount : 0) << '\n'
                     << "baseline_untracked_excitation_rate_s1="
                     << (baseEmission ? baseEmission->untrackedExcitationRate_s1 : 0.0)
                     << '\n'
                     << "anomaly_untracked_excitation_rate_s1="
                     << (changedEmission ? changedEmission->untrackedExcitationRate_s1 : 0.0)
                     << '\n'
                     << "cross_sections_root=" << (aurora ? iaaRoot.string() : "none") << '\n'
                     << "comparison_csv=" << comparisonPath.string() << '\n';
            if (!manifest) throw std::runtime_error("Comparison manifest write failed");
            std::cout << "Paired results: baseline " << baseline.precipitationCount()
                      << " entered, " << baseline.countOutcome(ParticleOutcome::Mirrored)
                      << " mirrored; anomaly " << changed.precipitationCount()
                      << " entered, " << changed.countOutcome(ParticleOutcome::Mirrored)
                      << " mirrored. Saved " << comparisonPath << '\n';
            return 0;
        }
        const ElectronSceneData scene = traceScene(
            baseField, initialConditions, threads, "Electron trajectories");
        std::cout << "Simulated " << scene.paths.size() << " markers; "
                  << scene.precipitationCount() << " precipitated, "
                  << scene.countOutcome(ParticleOutcome::Mirrored) << " mirrored, "
                  << scene.countOutcome(ParticleOutcome::Escaped) << " escaped, "
                  << scene.countOutcome(ParticleOutcome::Unresolved) << " unresolved"
                  << std::endl;
        if (!output.parent_path().empty()) {
            std::filesystem::create_directories(output.parent_path());
        }
        const auto outcomesPath = output.parent_path() /
            (output.stem().string() + "_particles.csv");
        saveParticleOutcomes(scene, outcomesPath);
        if (!framesDirectory.empty() &&
            std::filesystem::exists(framesDirectory) &&
            !std::filesystem::is_empty(framesDirectory)) {
            throw std::invalid_argument("--frames-dir must be empty");
        }
        VtkElectronScene visualizer(scene, width, height, !framesDirectory.empty());
        if (!emissionCsv.empty()) {
            const VolumeEmissionGrid emission = EmissionCsvReplay::load(emissionCsv);
            addGlow(visualizer, emission);
            std::cout << "Replaying saved emission grid " << emissionCsv << '\n';
        } else if (aurora && scene.precipitationCount() == 0) {
            std::cout << "No precipitating markers; skipping atmospheric emission" << std::endl;
        } else if (aurora) {
            const AuroraPipelineResult result = makeAurora(
                scene, iaaRoot, threads, seed, cascadeMarkers);
            addGlow(visualizer, result.emission);
            if (!output.parent_path().empty()) {
                std::filesystem::create_directories(output.parent_path());
            }
            const auto emissionPath = output.parent_path() /
                (output.stem().string() + "_emission.csv");
            result.emission.saveSparseCsv(emissionPath);
            const auto spectralDirectory = output.parent_path() /
                (output.stem().string() + "_spectral");
            saveSpectralViews(scene, result.emission, spectralDirectory);
            const auto manifestPath = output.parent_path() /
                (output.stem().string() + "_manifest.txt");
            std::ofstream manifest(manifestPath);
            if (!manifest) throw std::runtime_error("Cannot write scene manifest");
            manifest << std::setprecision(17)
                     << "seed=" << seed << '\n'
                     << "visual_markers=" << markerCount << '\n'
                     << "pitch_min_deg=" << minPitch_deg << '\n'
                     << "pitch_max_deg=" << maxPitch_deg << '\n'
                     << "precipitated_markers=" << scene.precipitationCount() << '\n'
                     << "mirrored_markers=" << scene.countOutcome(ParticleOutcome::Mirrored) << '\n'
                     << "escaped_markers=" << scene.countOutcome(ParticleOutcome::Escaped) << '\n'
                     << "unresolved_markers=" << scene.countOutcome(ParticleOutcome::Unresolved) << '\n'
                     << "particles_csv=" << outcomesPath.string() << '\n'
                     << "cascade_primaries=" << result.primaryCount << '\n'
                     << "threads=" << threads << '\n'
                     << "source_rate_electrons_s1=" << sourceRate_s1 << '\n'
                     << "atmosphere=NRLMSIS2.0 illustrative 2020-03-20 70N 0E\n"
                     << "cross_sections_root=" << iaaRoot.string() << '\n'
                     << "cross_section_tail=zero below 100eV endpoint; inverse-E above\n"
                     << "cross_section_model_max_eV=25000\n"
                     << "emission_grid_csv=" << emissionPath.string() << '\n'
                     << "spectral_directory=" << spectralDirectory.string() << '\n'
                     << "untracked_excitation_rate_s1="
                     << result.untrackedExcitationRate_s1 << '\n'
                     << "outside_grid_excitation_rate_s1="
                     << result.outsideGridRate_s1 << '\n'
                     << "energy_initial_eV_s1=" << result.initialEnergy_eV_s1 << '\n'
                     << "energy_residual_eV_s1=" << result.energyResidual_eV_s1 << '\n';
            if (!manifest) throw std::runtime_error("Scene manifest write failed");
            std::cout << "Atmospheric cascades: " << result.primaryCount
                      << " primaries, " << result.cascadeTracks << " tracks, "
                      << result.excitationEvents << " visible-line excitations; "
                      << result.untrackedExcitationRate_s1
                      << " untracked excitations/s; emission: "
                      << emissionPath << '\n';
        }
        if (!output.parent_path().empty()) {
            std::filesystem::create_directories(output.parent_path());
        }
        if (!framesDirectory.empty()) {
            std::filesystem::create_directories(framesDirectory);
            const double lastTime_s = ElectronAnimation::lastEventTime_s(scene);
            ProgressReporter frameProgress("Animation frames", frameCount, std::cerr);
            for (std::size_t frame = 0; frame < frameCount; ++frame) {
                const double time_s = ElectronAnimation::frameTime_s(
                    frame, frameCount, lastTime_s);
                visualizer.setAnimationFrame(scene, time_s);
                std::ostringstream name;
                name << "frame_" << std::setw(4) << std::setfill('0') << frame << ".png";
                visualizer.savePng(framesDirectory / name.str());
                frameProgress.advance();
            }
        }
        visualizer.savePng(output);
        std::cout << "Saved " << output << '\n';
        if (interactive) visualizer.interact();
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
