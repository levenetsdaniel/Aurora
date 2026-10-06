#ifndef AURORA_AURORAPIPELINE_HPP
#define AURORA_AURORAPIPELINE_HPP

#include "../atmosphere/ElectronCascade.hpp"
#include "../emission/VolumeEmissionGrid.hpp"
#include "../ensemble/PrecipitationCollector.hpp"
#include "../numerics/ParallelFor.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

struct AuroraPipelineConfig {
    ExcitationGridConfig grid;
    ElectronCascadeConfig cascade;
    std::size_t threads{1};
    std::size_t batchSize{64};
    std::uint64_t seed{1};
};

struct AuroraPipelineResult {
    VolumeEmissionGrid emission;
    std::size_t primaryCount{};
    std::size_t cascadeTracks{};
    std::size_t excitationEvents{};
    double untrackedExcitationRate_s1{};
    double outsideGridRate_s1{};
    double initialEnergy_eV_s1{};
    double energyResidual_eV_s1{};
};

class AuroraPipeline {
private:
    const AtmosphereModel& atmosphere;
    const CrossSectionDatabase& crossSections;
    const KineticsCatalog& states;
    const SpectralLineCatalog& lines;
    AuroraPipelineConfig config;

    static std::uint64_t markerSeed(std::uint64_t seed, std::uint64_t id) {
        // SplitMix64: the seed for one marker does not depend on thread scheduling.
        std::uint64_t value = seed + id + 0x9e3779b97f4a7c15ULL;
        value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
        value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
        return value ^ (value >> 31);
    }

    std::vector<ExcitedStateId> selectedStates() const {
        std::vector<ExcitedStateId> selected;
        for (const CollisionChannel& channel : crossSections.allChannels()) {
            if (channel.type != CollisionType::Excitation) continue;
            const ExcitedStateId state{channel.species, channel.targetState};
            try {
                const std::string& canonical = states.resolve(state).id;
                const bool visible = std::any_of(lines.lines().begin(), lines.lines().end(),
                    [&](const SpectralLine& line) {
                        return line.upperStateId == canonical;
                    });
                if (visible) selected.push_back(state);
            } catch (const std::out_of_range&) {
                // Other excitation channels remain collision-energy sinks.
            }
        }
        return selected;
    }

public:
    AuroraPipeline(const AtmosphereModel& atmosphere,
                   const CrossSectionDatabase& crossSections,
                   const KineticsCatalog& states,
                   const SpectralLineCatalog& lines,
                   AuroraPipelineConfig config)
        : atmosphere(atmosphere), crossSections(crossSections),
          states(states), lines(lines), config(config) {
        if (config.threads == 0 || config.batchSize == 0) {
            throw std::invalid_argument("AuroraPipeline: invalid parallel configuration");
        }
    }

    // Called once per completed primary, possibly from multiple worker threads.
    AuroraPipelineResult run(const std::vector<PrecipitatingElectron>& electrons,
                             const std::function<void()>& onPrimaryComplete = {}) const {
        const auto selected = selectedStates();
        if (selected.empty()) {
            throw std::invalid_argument("AuroraPipeline: no excitation channels match lines");
        }
        ExcitationModel excitation(crossSections, config.grid,
            ExcitationOutOfBounds::Ignore, selected);
        const ElectronTransportMC transport(atmosphere, crossSections);
        const ElectronCascade cascade(transport, config.cascade);
        AuroraPipelineResult summary{
            VolumeEmissionGrid(config.grid, lines.lines())};
        for (std::size_t begin = 0; begin < electrons.size(); begin += config.batchSize) {
            const std::size_t count = std::min(config.batchSize, electrons.size() - begin);
            std::vector<ElectronCascadeResult> results(count);
            parallelFor(count, config.threads, [&](std::size_t offset) {
                const PrecipitatingElectron& electron = electrons[begin + offset];
                const Vector3 velocity = Relativity::velocityFromMomentum(
                    electron.state.p, C::electronMass_kg);
                Random rng(markerSeed(config.seed, electron.id));
                results[offset] = cascade.run(electron.state.r, velocity,
                                              electron.energy_eV, rng);
                if (onPrimaryComplete) onPrimaryComplete();
            });
            // Deterministic merge, independent of scheduling and thread count.
            for (std::size_t offset = 0; offset < count; ++offset) {
                const PrecipitatingElectron& electron = electrons[begin + offset];
                const ElectronCascadeResult& result = results[offset];
                if (result.limitReached) {
                    throw std::runtime_error("AuroraPipeline: cascade limit for marker " +
                                             std::to_string(electron.id));
                }
                summary.initialEnergy_eV_s1 += result.ledger.initial_eV * electron.weight;
                summary.energyResidual_eV_s1 += result.ledger.residual_eV() * electron.weight;
                summary.cascadeTracks += result.tracks.size();
                summary.excitationEvents += excitation.accumulate(result, electron.weight);
                ++summary.primaryCount;
            }
        }
        summary.untrackedExcitationRate_s1 = excitation.untrackedEventRate_s1();
        summary.outsideGridRate_s1 = excitation.ignoredEventRate_s1();
        const ExcitedStateKinetics kinetics(states, atmosphere);
        summary.emission = VolumeEmissionGrid::fromKinetics(
            excitation, kinetics, states, lines);
        return summary;
    }
};

#endif // AURORA_AURORAPIPELINE_HPP
