#include "ElectronSceneData.hpp"
#include "VtkElectronScene.hpp"
#include "AuroraPipeline.hpp"
#include "../io/AtmosphereCsvLoader.hpp"
#include "../io/IaaCrossSectionLoader.hpp"
#include "../io/KineticsCsvLoader.hpp"
#include "../io/SpectralLineCsvLoader.hpp"
#include "../emission/OpticalRenderer.hpp"
#include "../numerics/Random.hpp"
#include "../numerics/ParallelFor.hpp"
#include "../physics/DipoleField.hpp"
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
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#ifndef AURORA_SOURCE_DIR
#error AURORA_SOURCE_DIR is required
#endif

namespace {

InitialCondition sampleElectron(const DipoleField& field, Random& random,
                                 std::uint64_t id, double markerRate_s1,
                                 bool auroralArc) {
    const double longitude = auroralArc ?
        -0.7 + 1.4 * random.uniform01() :
        2.0 * std::numbers::pi * random.uniform01();
    const double colatitude = auroralArc ?
        0.27 + 0.05 * random.uniform01() :
        0.20 + 0.17 * random.uniform01();
    const double radius_m = C::earthRadius_m + 450000.0;
    const Vector3 position_m(
        radius_m * std::sin(colatitude) * std::cos(longitude),
        radius_m * std::sin(colatitude) * std::sin(longitude),
        radius_m * std::cos(colatitude));
    const Vector3 magnetic = field.evaluate(position_m);
    const Vector3 along = magnetic * (-1.0 / magnetic.norm());
    const Vector3 reference(0.0, 0.0, 1.0);
    const Vector3 across = (reference - along * reference.dot(along)).normalize();
    const Vector3 around = along.cross(across);
    const double pitch = 0.04 + 0.18 * random.uniform01();
    const double phase = 2.0 * std::numbers::pi * random.uniform01();
    const Vector3 direction = along * std::cos(pitch) +
        (across * std::cos(phase) + around * std::sin(phase)) * std::sin(pitch);
    const double energy_eV = 5000.0 + 20000.0 * random.uniform01();
    const double momentum = Relativity::momentumMagnitudeFromKineticEnergy_J(
        energy_eV * C::elementaryCharge_C, C::electronMass_kg);
    return {id, {position_m, direction * momentum, 0.0}, energy_eV,
            pitch, phase, markerRate_s1};
}

ElectronSceneData makeScene(std::size_t markerCount, std::uint64_t seed,
                            std::size_t threads, double sourceRate_s1,
                            bool auroralArc) {
    constexpr double atmosphereHeight_m = 150000.0;
    DipoleField field;
    const SimulationConfig config{2e-8, atmosphereHeight_m,
        C::earthRadius_m + 800000.0, 0.015, 750000};
    const SimulationEngine engine(field, -C::elementaryCharge_C,
                                  C::electronMass_kg, config);
    const ElectronSceneSampler sampler(engine, field);
    Random random(seed);
    ElectronSceneData scene;
    scene.atmosphereHeight_m = atmosphereHeight_m;
    scene.paths.reserve(markerCount);
    const double markerRate_s1 = sourceRate_s1 / static_cast<double>(markerCount);
    std::vector<InitialCondition> initialConditions;
    initialConditions.reserve(markerCount);
    for (std::size_t index = 0; index < markerCount; ++index) {
        initialConditions.push_back(sampleElectron(
            field, random, index + 1, markerRate_s1, auroralArc));
    }
    scene.paths.resize(markerCount);
    parallelFor(markerCount, threads, [&](std::size_t index) {
        scene.paths[index] = sampler.traceElectron(initialConditions[index], 250);
    });
    for (std::size_t index = 0; index < markerCount; ++index) {
        if (index % std::max<std::size_t>(1, markerCount / 12) == 0) {
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
                                std::size_t cascadeMarkers) {
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
    if (cascadeMarkers > 0 && cascadeMarkers < electrons.size()) {
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
    return pipeline.run(electrons);
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

} // namespace

int main(int argc, char** argv) {
    try {
        std::filesystem::path output = "electron_flow.png";
        std::size_t markerCount = 48;
        std::size_t cascadeMarkers = 0;
        std::uint64_t seed = 20260928;
        std::size_t threads = std::max(1u, std::thread::hardware_concurrency());
        double sourceRate_s1 = 1.0e20;
        std::filesystem::path iaaRoot;
        bool aurora = false;
        bool interactive = false;
        for (int index = 1; index < argc; ++index) {
            const std::string argument = argv[index];
            if (argument == "--output" && index + 1 < argc) {
                output = argv[++index];
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
            } else if (argument == "--iaa-root" && index + 1 < argc) {
                iaaRoot = argv[++index];
            } else if (argument == "--aurora") {
                aurora = true;
            } else if (argument == "--interactive") {
                interactive = true;
            } else {
                throw std::invalid_argument(
                    "Usage: aurora_vtk_visualizer [--output image.png] "
                    "[--markers N] [--cascade-markers N] [--seed N] "
                    "[--threads N] [--source-rate N] "
                    "[--aurora --iaa-root DIR] [--interactive]");
            }
        }
        if (markerCount == 0 || markerCount > 5000 ||
            threads == 0 || threads > 64 ||
            !std::isfinite(sourceRate_s1) || sourceRate_s1 <= 0.0 ||
            output.empty() || output.extension() != ".png" ||
            (aurora && iaaRoot.empty())) {
            throw std::invalid_argument("Invalid marker count, rate, thread count, IAA root or output path");
        }
        const ElectronSceneData scene = makeScene(
            markerCount, seed, threads, sourceRate_s1, aurora);
        std::cout << "Simulated " << scene.paths.size() << " markers; "
                  << scene.precipitationCount() << " reached the atmosphere\n";
        VtkElectronScene visualizer(scene);
        if (aurora) {
            const AuroraPipelineResult result = makeAurora(
                scene, iaaRoot, threads, seed, cascadeMarkers);
            double greenPeak = 0.0;
            double redPeak = 0.0;
            const auto& grid = result.emission.configuration();
            for (std::size_t z = 0; z < grid.cellsZ; ++z) {
                for (std::size_t y = 0; y < grid.cellsY; ++y) {
                    for (std::size_t x = 0; x < grid.cellsX; ++x) {
                        greenPeak = std::max(greenPeak,
                            result.emission.emissivity("O_5577", x, y, z));
                        redPeak = std::max(redPeak,
                            result.emission.emissivity("O_6300", x, y, z));
                    }
                }
            }
            const double referencePeak = std::max(greenPeak, redPeak);
            const double greenOpacity = referencePeak > 0.0 ?
                std::sqrt(greenPeak / referencePeak) : 0.0;
            const double redOpacity = referencePeak > 0.0 ?
                std::sqrt(redPeak / referencePeak) : 0.0;
            visualizer.setCameraScale(0.31);
            visualizer.addEmission(result.emission, "O_6300", 1.0, 0.16, 0.10,
                                   redOpacity);
            visualizer.addEmission(result.emission, "O_5577", 0.14, 1.0, 0.27,
                                   greenOpacity);
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
        visualizer.savePng(output);
        std::cout << "Saved " << output << '\n';
        if (interactive) visualizer.interact();
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
