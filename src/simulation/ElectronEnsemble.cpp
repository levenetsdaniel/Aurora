#include "aurora/simulation/ElectronEnsemble.hpp"
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>

namespace aurora {
std::vector<ElectronResult> run_ensemble(
    const EnsembleParameters& parameters, const std::vector<ElectronSample>& samples) {
    if (parameters.particle_count == 0)
        throw std::invalid_argument("Ensemble must contain at least one electron");
    if (parameters.particle_count != samples.size())
        throw std::invalid_argument("particle_count must match the number of samples");

    std::unordered_set<std::size_t> identifiers;
    identifiers.reserve(samples.size());
    std::vector<PreparedExperiment> prepared;
    prepared.reserve(samples.size());

    for (const auto& sample : samples) {
        if (!identifiers.insert(sample.id).second)
            throw std::invalid_argument("Duplicate electron id: " + std::to_string(sample.id));

        auto particle = parameters.particle;
        particle.energy_keV = sample.energy_keV;
        particle.pitch_deg = sample.pitch_deg;
        particle.gyro_phase_deg = sample.gyro_phase_deg;

        prepared.push_back(prepare_experiment(particle));
    }

    std::vector<ElectronResult> results;
    results.reserve(samples.size());

    for (std::size_t index = 0; index < samples.size(); ++index) {
        auto simulation = run_atmosphere_experiment(prepared[index]);
        results.push_back({samples[index], std::move(simulation)});
    }

    return results;
}
}
