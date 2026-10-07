#include "PrecipitationStatistics.hpp"
#include "../physics/LocalizedDipoleField.hpp"

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
#include <vector>

namespace {

constexpr double colatitude_rad = 0.32;
constexpr double sourceHeight_m = 450000.0;
constexpr double boundaryHeight_m = 150000.0;
constexpr double degrees = 180.0 / std::numbers::pi;

struct Case {
    const char* name;
    double depth_m;
    double surfaceFraction;
};

constexpr Case cases[]{
    {"zero", 500000.0, 0.0},
    {"weak_broad", 500000.0, 0.01},
    {"strong_broad", 500000.0, 0.30},
    {"strong_compact", 250000.0, 0.30},
};

const char* outcomeName(ParticleOutcome outcome) {
    switch (outcome) {
        case ParticleOutcome::Precipitated: return "Precipitated";
        case ParticleOutcome::Mirrored: return "Mirrored";
        case ParticleOutcome::Escaped: return "Escaped";
        case ParticleOutcome::Unresolved: return "Unresolved";
    }
    throw std::logic_error("Unknown outcome");
}

std::ofstream outputFile(const std::filesystem::path& directory,
                         const char* name) {
    std::ofstream output(directory / name);
    if (!output) throw std::runtime_error(std::string("Cannot write ") + name);
    output << std::setprecision(17);
    return output;
}

Vector3 radialPoint(double altitude_m, double offset_rad) {
    const double angle = colatitude_rad + offset_rad;
    return {(C::earthRadius_m + altitude_m) * std::sin(angle), 0.0,
            (C::earthRadius_m + altitude_m) * std::cos(angle)};
}

double surfacePeak(const Case& sample, const DipoleField& background) {
    return sample.surfaceFraction *
           background.evaluate(radialPoint(0.0, 0.0)).norm();
}

struct Bracket {
    double low_deg;
    double high_deg;
};

// All cases use the same BASE-FIELD pitch axis and gyro phase. The anomaly
// must never silently redefine the source distribution.
Bracket scan(const MagneticField& field, const Vector3& position,
             const Vector3& along, const Vector3& across, double momentum,
             double dt_s) {
    const SimulationEngine engine(field, -C::elementaryCharge_C, C::electronMass_kg,
        {dt_s, boundaryHeight_m, C::earthRadius_m + 800000.0, 0.05,
         static_cast<std::uint64_t>(std::ceil(0.05 / dt_s)) + 1});
    const auto outcome = [&](double pitch_deg) {
        const double pitch = pitch_deg / degrees;
        const ParticleState initial{position,
            (along * std::cos(pitch) + across * std::sin(pitch)) * momentum, 0.0};
        return engine.run(initial).outcome;
    };
    double low = 50.0;
    double high = 85.0;
    if (outcome(low) != ParticleOutcome::Precipitated ||
        outcome(high) != ParticleOutcome::Mirrored)
        throw std::runtime_error("Anomaly scan did not bracket a loss cone");
    while (high - low > 0.05) {
        const double middle = (low + high) / 2.0;
        const auto result = outcome(middle);
        if (result == ParticleOutcome::Precipitated) low = middle;
        else if (result == ParticleOutcome::Mirrored) high = middle;
        else throw std::runtime_error("Anomaly scan reached Escaped/Unresolved");
    }
    return {low, high};
}

struct ParticleComparison {
    SimulationResult result;
    double maxRelativeAnomaly{};
};

} // namespace

int main(int argc, char* argv[]) {
    try {
        if (argc != 2 && argc != 4) {
            std::cerr << "Usage: aurora_anomaly_experiment OUTPUT_DIRECTORY [--threads N]\n";
            return 2;
        }
        std::size_t threads = std::min<std::size_t>(8,
            std::max(1u, std::thread::hardware_concurrency()));
        if (argc == 4) {
            if (std::string(argv[2]) != "--threads")
                throw std::invalid_argument("Expected --threads N");
            threads = std::stoull(argv[3]);
            if (threads == 0 || threads > 64)
                throw std::invalid_argument("Thread count must be in [1,64]");
        }
        const std::filesystem::path directory(argv[1]);
        std::filesystem::create_directories(directory);
        auto fieldCsv = outputFile(directory, "anomaly_field.csv");
        auto scanCsv = outputFile(directory, "anomaly_scan.csv");
        auto summaryCsv = outputFile(directory, "anomaly_summary.csv");
        auto particleCsv = outputFile(directory, "anomaly_particles.csv");
        fieldCsv << "case,depth_m,surface_fraction,surface_peak_T,altitude_m,offset_deg,"
                    "background_B_T,anomaly_B_T,relative_anomaly,anomaly_radial_T\n";
        scanCsv << "case,depth_m,surface_fraction,dt_s,baseline_low_deg,"
                   "baseline_high_deg,anomaly_low_deg,anomaly_high_deg,"
                   "delta_midpoint_deg,delta_midpoint_fine_deg\n";
        summaryCsv << "case,depth_m,surface_fraction,seed,N,dt_s,baseline_precipitated,"
                      "baseline_mirrored,baseline_escaped,baseline_unresolved,"
                      "anomaly_precipitated,anomaly_mirrored,anomaly_escaped,"
                      "anomaly_unresolved,precip_to_mirror,mirror_to_precip,"
                      "other_transitions,mean_entry_shift_m,max_relative_anomaly\n";
        particleCsv << "case,depth_m,surface_fraction,seed,id,energy_eV,pitch_deg,"
                       "phase_deg,initial_time_s,initial_x_m,initial_y_m,initial_z_m,"
                       "initial_px_kg_m_s,initial_py_kg_m_s,initial_pz_kg_m_s,"
                       "sample_weight,baseline_outcome,anomaly_outcome,"
                       "baseline_event_x_m,baseline_event_y_m,baseline_event_z_m,"
                       "anomaly_event_x_m,anomaly_event_y_m,anomaly_event_z_m,"
                       "entry_shift_m,max_relative_anomaly\n";

        const DipoleField background;
        const Vector3 position = radialPoint(sourceHeight_m, 0.0);
        const Vector3 baseB = background.evaluate(position);
        const Vector3 along = baseB * (-1.0 / baseB.norm());
        const Vector3 across(0.0, 1.0, 0.0);
        const double momentum = TrajectoryExperiments::momentumMagnitude();
        const double dt_s = TrajectoryExperiments::gyroPeriod(baseB.norm()) / 50.0;
        const Bracket baselineBracket = scan(background, position, along, across,
                                              momentum, dt_s);
        const double critical = (baselineBracket.low_deg + baselineBracket.high_deg) / 2.0;
        const std::uint64_t seeds[]{20261006, 20261007};
        std::vector<PrecipitationStatistics::Result> baselineRuns;
        for (const auto seed : seeds) {
            const PrecipitationStatistics::Config config{
                colatitude_rad, PrecipitationStatistics::PitchDistribution::NarrowBin,
                critical, 1.2, 128, seed, threads, 50};
            baselineRuns.push_back(PrecipitationStatistics::run(config));
        }

        for (const auto& sample : cases) {
            std::cout << "Computing anomaly case " << sample.name << std::endl;
            const double peak = surfacePeak(sample, background);
            const LocalizedDipoleField field(colatitude_rad, sample.depth_m, peak);
            for (double altitude : {0.0, 150000.0, 450000.0}) {
                for (int offset = -12; offset <= 12; ++offset) {
                    const Vector3 point = radialPoint(altitude, offset * 0.5 / degrees);
                    const Vector3 base = background.evaluate(point);
                    const Vector3 delta = field.anomaly(point);
                    const Vector3 radial = point / point.norm();
                    fieldCsv << sample.name << ',' << sample.depth_m << ','
                             << sample.surfaceFraction << ',' << peak << ','
                             << altitude << ',' << offset * 0.5 << ','
                             << base.norm() << ',' << delta.norm() << ','
                             << delta.norm() / base.norm() << ','
                             << delta.dot(radial) << '\n';
                }
            }
            const Bracket anomalyBracket = scan(field, position, along, across,
                                                momentum, dt_s);
            const Bracket fineBaseline = scan(background, position, along, across,
                                              momentum, dt_s / 2.0);
            const Bracket fineAnomaly = scan(field, position, along, across,
                                             momentum, dt_s / 2.0);
            scanCsv << sample.name << ',' << sample.depth_m << ','
                    << sample.surfaceFraction << ',' << dt_s << ','
                    << baselineBracket.low_deg << ',' << baselineBracket.high_deg << ','
                    << anomalyBracket.low_deg << ',' << anomalyBracket.high_deg << ','
                    << (anomalyBracket.low_deg + anomalyBracket.high_deg -
                        baselineBracket.low_deg - baselineBracket.high_deg) / 2.0
                    << ',' << (fineAnomaly.low_deg + fineAnomaly.high_deg -
                               fineBaseline.low_deg - fineBaseline.high_deg) / 2.0
                    << '\n';

            for (std::size_t seedIndex = 0; seedIndex < 2; ++seedIndex) {
                const auto& baseline = baselineRuns[seedIndex];
                const SimulationEngine engine(field, -C::elementaryCharge_C,
                    C::electronMass_kg,
                    {dt_s, boundaryHeight_m, C::earthRadius_m + 800000.0, 0.05,
                     static_cast<std::uint64_t>(std::ceil(0.05 / dt_s)) + 1});
                std::vector<ParticleComparison> compared(baseline.initialConditions.size());
                parallelFor(compared.size(), threads, [&](std::size_t index) {
                    std::size_t observation = 0;
                    double maximum = 0.0;
                    const auto observe = [&](const ParticleState& state) {
                        if (observation++ % 50 != 0) return;
                        const Vector3 base = background.evaluate(state.r);
                        maximum = std::max(maximum,
                            field.anomaly(state.r).norm() / base.norm());
                    };
                    compared[index].result = engine.run(
                        baseline.initialConditions[index].state, observe);
                    const Vector3 end = compared[index].result.eventState.r;
                    maximum = std::max(maximum,
                        field.anomaly(end).norm() / background.evaluate(end).norm());
                    compared[index].maxRelativeAnomaly = maximum;
                });
                const auto anomalyResults = [&] {
                    std::vector<SimulationResult> results;
                    results.reserve(compared.size());
                    for (const auto& item : compared) results.push_back(item.result);
                    return results;
                }();
                const auto aCounts = PrecipitationStatistics::count(
                    anomalyResults, anomalyResults.size());
                std::size_t precipToMirror = 0;
                std::size_t mirrorToPrecip = 0;
                std::size_t otherTransitions = 0;
                std::size_t commonEntries = 0;
                double entryShiftSum = 0.0;
                double maximum = 0.0;
                for (std::size_t index = 0; index < compared.size(); ++index) {
                    const auto& initial = baseline.initialConditions[index];
                    const auto& first = baseline.trajectories[index];
                    const auto& second = compared[index].result;
                    if (first.outcome == ParticleOutcome::Precipitated &&
                        second.outcome == ParticleOutcome::Mirrored) ++precipToMirror;
                    else if (first.outcome == ParticleOutcome::Mirrored &&
                             second.outcome == ParticleOutcome::Precipitated) ++mirrorToPrecip;
                    else if (first.outcome != second.outcome) ++otherTransitions;
                    double entryShift = 0.0;
                    if (first.outcome == ParticleOutcome::Precipitated &&
                        second.outcome == ParticleOutcome::Precipitated) {
                        entryShift = (second.eventState.r - first.eventState.r).norm();
                        entryShiftSum += entryShift;
                        ++commonEntries;
                    }
                    maximum = std::max(maximum, compared[index].maxRelativeAnomaly);
                    particleCsv << sample.name << ',' << sample.depth_m << ','
                                << sample.surfaceFraction << ',' << seeds[seedIndex] << ','
                                << initial.id << ',' << initial.energy_eV << ','
                                << initial.pitch_rad * degrees << ','
                                << initial.gyroPhase_rad * degrees << ','
                                << initial.state.time_s << ','
                                << initial.state.r.x << ',' << initial.state.r.y << ','
                                << initial.state.r.z << ',' << initial.state.p.x << ','
                                << initial.state.p.y << ',' << initial.state.p.z << ','
                                << initial.weight << ',' << outcomeName(first.outcome) << ','
                                << outcomeName(second.outcome) << ','
                                << first.eventState.r.x << ',' << first.eventState.r.y << ','
                                << first.eventState.r.z << ',' << second.eventState.r.x << ','
                                << second.eventState.r.y << ',' << second.eventState.r.z << ','
                                << entryShift << ',' << compared[index].maxRelativeAnomaly
                                << '\n';
                    if (sample.surfaceFraction == 0.0 &&
                        (first.outcome != second.outcome ||
                         (first.eventState.r - second.eventState.r).norm() != 0.0))
                        throw std::runtime_error("Zero anomaly changed a paired trajectory");
                }
                summaryCsv << sample.name << ',' << sample.depth_m << ','
                           << sample.surfaceFraction << ',' << seeds[seedIndex] << ','
                           << compared.size() << ',' << dt_s << ','
                           << baseline.counts.precipitated << ','
                           << baseline.counts.mirrored << ',' << baseline.counts.escaped
                           << ',' << baseline.counts.unresolved << ','
                           << aCounts.precipitated << ',' << aCounts.mirrored << ','
                           << aCounts.escaped << ',' << aCounts.unresolved << ','
                           << precipToMirror << ',' << mirrorToPrecip << ','
                           << otherTransitions << ','
                           << (commonEntries ? entryShiftSum / commonEntries : 0.0)
                           << ',' << maximum << '\n';
            }
            fieldCsv.flush();
            scanCsv.flush();
            summaryCsv.flush();
            particleCsv.flush();
        }
        if (!fieldCsv || !scanCsv || !summaryCsv || !particleCsv)
            throw std::runtime_error("Anomaly CSV write failed");
        std::cout << "Saved anomaly experiment CSV files to " << directory << '\n';
    } catch (const std::exception& error) {
        std::cerr << "Anomaly experiment failed: " << error.what() << '\n';
        return 1;
    }
}
