#include "TrajectoryExperiments.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

std::ofstream outputFile(const std::filesystem::path& directory,
                         const std::string& name) {
    std::ofstream output(directory / name);
    if (!output) throw std::runtime_error("cannot create " + (directory / name).string());
    output << std::setprecision(17);
    return output;
}

void writeTrace(std::ofstream& output, const std::string& field,
                const std::vector<TrajectoryExperiments::DiagnosticRow>& rows) {
    for (const auto& row : rows) {
        output << field << ',' << row.time_s << ',' << row.energy_eV << ','
               << row.relativeEnergyError << ',' << row.magneticMagnitude_T << ','
               << row.parallelMomentum_kg_m_s << ',' << row.magneticMoment_J_T << '\n';
    }
}

const char* outcomeName(ParticleOutcome outcome) {
    switch (outcome) {
        case ParticleOutcome::Precipitated: return "Precipitated";
        case ParticleOutcome::Mirrored: return "Mirrored";
        case ParticleOutcome::Escaped: return "Escaped";
        case ParticleOutcome::Unresolved: return "Unresolved";
    }
    throw std::runtime_error("unknown particle outcome");
}

} // namespace

int main(int argc, char* argv[]) {
    try {
        if (argc != 2) {
            std::cerr << "Usage: aurora_trajectory_experiments OUTPUT_DIRECTORY\n";
            return 2;
        }
        const std::filesystem::path directory(argv[1]);
        std::filesystem::create_directories(directory);

        auto orbit = outputFile(directory, "uniform_orbit.csv");
        orbit << "steps_per_period,dt_s,exact_period_s,measured_period_s,exact_radius_m,"
                 "max_radius_error_m,position_error_m,max_relative_energy_error\n";
        auto event = outputFile(directory, "atmosphere_event.csv");
        event << "steps_per_period,dt_s,exact_time_s,simulated_time_s,time_error_s,"
                 "sphere_residual_m,steps\n";
        for (int resolution : {50, 100, 200, 400}) {
            const auto o = TrajectoryExperiments::uniformOrbit(resolution);
            orbit << o.stepsPerPeriod << ',' << o.dt_s << ',' << o.exactPeriod_s << ','
                  << o.measuredPeriod_s << ',' << o.exactRadius_m << ','
                  << o.radiusError_m << ',' << o.positionError_m << ','
                  << o.maximumRelativeEnergyError << '\n';
            const auto e = TrajectoryExperiments::atmosphericCrossing(resolution);
            event << e.stepsPerPeriod << ',' << e.dt_s << ',' << e.exactTime_s << ','
                  << e.simulatedTime_s << ',' << e.timeError_s << ','
                  << e.sphereResidual_m << ',' << e.steps << '\n';
        }

        const double momentum = TrajectoryExperiments::momentumMagnitude();
        const double uniformDt = TrajectoryExperiments::gyroPeriod(
            TrajectoryExperiments::field_T) / 100.0;
        const UniformField uniform(Vector3(0.0, 0.0, TrajectoryExperiments::field_T));
        const ParticleState uniformInitial{Vector3(C::earthRadius_m + 450000.0, 0.0, 0.0),
            Vector3(0.8 * momentum, 0.0, 0.6 * momentum), 0.0};
        const DipoleField dipole;
        const Vector3 dipolePosition(C::earthRadius_m + 450000.0, 0.0, 0.0);
        const Vector3 dipoleField = dipole.evaluate(dipolePosition);
        const Vector3 dipoleAlong = dipoleField / dipoleField.norm();
        const ParticleState dipoleInitial{dipolePosition,
            (Vector3(1.0, 0.0, 0.0) * 0.8 + dipoleAlong * 0.6) * momentum, 0.0};
        const double dipoleDt = TrajectoryExperiments::gyroPeriod(
            dipole.evaluate(dipolePosition).norm()) / 100.0;
        auto energy = outputFile(directory, "energy_trace.csv");
        energy << "field,time_s,energy_eV,relative_energy_error,B_T,"
                  "p_parallel_kg_m_s,mu_J_T\n";
        writeTrace(energy, "uniform", TrajectoryExperiments::energyTrace(
            uniform, uniformInitial, uniformDt, 4000, 20));
        writeTrace(energy, "dipole", TrajectoryExperiments::energyTrace(
            dipole, dipoleInitial, dipoleDt, 4000, 20));

        const auto mirror = TrajectoryExperiments::mirrorTrajectory();
        auto mirrorSummary = outputFile(directory, "mirror_summary.csv");
        mirrorSummary << "outcome,B_initial_T,B_event_T,B_m_predicted_T,"
                         "relative_energy_error,steps\n";
        mirrorSummary << outcomeName(mirror.outcome) << ','
                      << mirror.initialMagneticMagnitude_T << ','
                      << mirror.mirrorMagneticMagnitude_T << ','
                      << mirror.predictedMirrorMagneticMagnitude_T << ','
                      << mirror.relativeEnergyError << ',' << mirror.steps << '\n';
        auto mirrorTrace = outputFile(directory, "mirror_trace.csv");
        mirrorTrace << "field,time_s,energy_eV,relative_energy_error,B_T,"
                       "p_parallel_kg_m_s,mu_J_T\n";
        writeTrace(mirrorTrace, "dipole", mirror.trace);

        std::cout << "Wrote experiment CSV files to " << directory << '\n';
        std::cout << "Mirror outcome: " << outcomeName(mirror.outcome)
                  << ", steps: " << mirror.steps << '\n';
        return mirror.outcome == ParticleOutcome::Mirrored ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "Trajectory experiment failed: " << error.what() << '\n';
        return 1;
    }
}
