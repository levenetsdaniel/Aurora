#include "aurora/io/ExperimentFiles.hpp"
#include "aurora/physics/Relativity.hpp"
#include <array>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numbers>
#include <sstream>
#include <vector>

namespace {
using namespace aurora;

struct Case {
    ExperimentParameters parameters;
    SavedExperiment saved;
    std::string kind;
    double offset_deg;
};

struct Bracket {
    double dt, lower, upper;
    int hemisphere;
    std::size_t lower_id, upper_id;
};

std::filesystem::path session_directory(const std::filesystem::path& root) {
    std::filesystem::create_directories(root);
    const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now().time_since_epoch())
                        .count();
    for (int i = 0; i < 10000; ++i) {
        const auto path = std::filesystem::absolute(
            root / ("precipitation-" + std::to_string(ms) + "-" + std::to_string(i)));
        if (std::filesystem::create_directory(path))
            return path;
    }
    throw std::runtime_error("Cannot create experiment session");
}
}

int main(int argc, char** argv) {
    try {
        std::filesystem::path root = "results/precipitation";
        if (argc == 3 && std::string(argv[1]) == "--output")
            root = argv[2];
        else if (argc != 1)
            throw std::invalid_argument(
                "Usage: aurora_atmospheric_precipitation [--output directory]");
        const auto directory = session_directory(root);
        std::ofstream csv(directory / "cases.csv");
        csv.exceptions(std::ios::failbit | std::ios::badbit);
        csv << std::setprecision(17)
            << "case_id,kind,energy_keV,pitch_deg,gyro_phase_deg,offset_from_cone_deg,dt_s,prediction,outcome,agreement,"
               "time_s,steps,altitude_m,x_m,y_m,z_m,max_energy_error,max_mu_error,wall_time_s,run_directory\n";
        std::vector<Case> cases;
        std::map<std::string, std::size_t> cache;
        ExperimentParameters base;
        base.max_time_s = 0.5;
        base.max_steps = 100000000;
        base.sample_every = 5000;
        const auto e = prepare_experiment(base);
        const double critical = e.loss_cone.critical_angle_rad * 180 / std::numbers::pi;
        bool passed = true;
        auto run_case = [&](ExperimentParameters p, const std::string& kind, double offset) {
            std::ostringstream key;
            write_experiment_parameters(key, p);
            if (const auto found = cache.find(key.str()); found != cache.end())
                return found->second;
            auto saved = save_experiment(p, directory / "cases");
            const auto prepared = prepare_experiment(p);
            const auto& result = saved.result;
            const auto& state = result.trajectory.final_state;
            const std::size_t id = cases.size();
            csv << id << ',' << kind << ',' << p.energy_keV << ',' << p.pitch_deg << ','
                << p.gyro_phase_deg << ',' << offset << ',' << p.dt_s << ','
                << (prepared.prediction ? prediction_name(*prepared.prediction) : "not_applicable")
                << ',' << outcome_name(result.outcome) << ',' << agreement_name(result.agreement)
                << ',' << state.time_s << ',' << result.trajectory.completed_steps << ','
                << norm(state.position_m) - constants::earth_radius_m << ',' << state.position_m.x
                << ',' << state.position_m.y << ',' << state.position_m.z << ','
                << result.trajectory.max_relative_energy_error << ','
                << result.trajectory.max_relative_magnetic_moment_error << ',' << saved.wall_time_s
                << ',' << std::filesystem::relative(saved.directory, directory).generic_string()
                << '\n';
            csv.flush();
            if (kind != "time_limit" && kind != "step_limit") {
                passed = passed && is_resolved(result.outcome) &&
                    result.outcome != ExperimentOutcome::Escaped &&
                    result.trajectory.max_relative_energy_error < 1e-9;
                if (std::abs(offset) >= 0.099 && prepared.prediction)
                    passed = passed && result.agreement == PredictionAgreement::Match;
            }
            std::cout << "Case " << id << ": " << p.energy_keV << " keV, " << p.pitch_deg
                      << " deg, dt=" << p.dt_s << " -> " << outcome_name(result.outcome) << ' '
                      << result.trajectory.error_message << '\n'
                      << std::flush;
            cases.push_back({p, std::move(saved), kind, offset});
            cache.emplace(key.str(), id);
            return id;
        };
        auto parameters = [&](double energy, double local_angle, double dt, int hemisphere) {
            auto p = base;
            p.energy_keV = energy;
            p.dt_s = dt;
            p.pitch_deg = hemisphere > 0 ? local_angle : 180 - local_angle;
            return p;
        };
        for (double energy : {10., 100., 1000.})
            for (int hemisphere : {1, -1}) {
                for (double local_angle : {5., critical - 0.1, critical, critical + 0.1, 45.})
                    run_case(parameters(energy, local_angle, 4e-8, hemisphere), "scan",
                        local_angle - critical);
            }
        for (int hemisphere : {1, -1})
            for (double dt : {2e-8, 1e-8})
                for (double offset : {-0.1, 0., 0.1})
                    run_case(parameters(100, critical + offset, dt, hemisphere), "step_refinement",
                        offset);

        std::vector<Bracket> brackets;
        for (int hemisphere : {1, -1})
            for (double dt : {4e-8, 2e-8, 1e-8}) {
                double lo = critical - 0.1, hi = critical + 0.1;
                auto lo_id =
                    run_case(parameters(100, lo, dt, hemisphere), "boundary_search", lo - critical);
                auto hi_id =
                    run_case(parameters(100, hi, dt, hemisphere), "boundary_search", hi - critical);
                bool valid =
                    cases[lo_id].saved.result.outcome == ExperimentOutcome::AtmosphereReached &&
                    cases[hi_id].saved.result.outcome == ExperimentOutcome::Reflected;
                for (int iteration = 0; iteration < 5 && valid; ++iteration) {
                    const double mid = lo + (hi - lo) * 0.5;
                    const auto id = run_case(
                        parameters(100, mid, dt, hemisphere), "boundary_search", mid - critical);
                    const auto outcome = cases[id].saved.result.outcome;
                    if (outcome == ExperimentOutcome::AtmosphereReached) {
                        lo = mid;
                        lo_id = id;
                    } else if (outcome == ExperimentOutcome::Reflected) {
                        hi = mid;
                        hi_id = id;
                    } else
                        valid = false;
                }
                passed = passed && valid;
                brackets.push_back({dt, lo, hi, hemisphere, lo_id, hi_id});
            }
        for (int hemisphere : {1, -1})
            for (double phase : {90., 180., 270.}) {
                auto p = parameters(1000, critical, 2e-8, hemisphere);
                p.gyro_phase_deg = phase;
                run_case(p, "gyro_phase", 0);
            }
        auto limited = base;
        limited.max_time_s = 1e-6;
        const auto time_id = run_case(limited, "time_limit", limited.pitch_deg - critical);
        limited = base;
        limited.max_steps = 2;
        const auto steps_id = run_case(limited, "step_limit", limited.pitch_deg - critical);
        passed = passed && cases[time_id].saved.result.outcome == ExperimentOutcome::TimeLimit &&
            cases[steps_id].saved.result.outcome == ExperimentOutcome::StepLimit;
        csv.close();

        std::array<std::size_t, 6> counts{};
        for (const auto& c : cases)
            ++counts[static_cast<std::size_t>(c.saved.result.outcome)];

        std::ofstream summary(directory / "summary.csv");
        summary.exceptions(std::ios::failbit | std::ios::badbit);
        summary
            << "passed,cases,atmosphere_reached,reflected,escaped,time_limit,step_limit,numerical_error\n"
            << (passed ? "true" : "false") << ',' << cases.size();
        for (auto count : counts)
            summary << ',' << count;
        summary << '\n';
        summary.close();

        std::ofstream boundaries(directory / "boundaries.csv");
        boundaries.exceptions(std::ios::failbit | std::ios::badbit);
        boundaries << std::setprecision(17)
                   << "direction,dt_s,entry_angle_deg,reflection_angle_deg,width_deg,"
                      "analytical_deg,entry_case_id,reflection_case_id,valid\n";
        for (const auto& b : brackets) {
            const bool valid =
                cases[b.lower_id].saved.result.outcome == ExperimentOutcome::AtmosphereReached &&
                cases[b.upper_id].saved.result.outcome == ExperimentOutcome::Reflected &&
                b.upper - b.lower <= 0.0062500001;
            boundaries << b.hemisphere << ',' << b.dt << ',' << b.lower << ',' << b.upper << ','
                       << b.upper - b.lower << ',' << critical << ',' << b.lower_id << ','
                       << b.upper_id << ',' << (valid ? "true" : "false") << '\n';
        }
        boundaries.close();

        std::cout << (passed ? "PASS: " : "FAIL: ") << directory.string() << '\n';
        return passed ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "precipitation: " << error.what() << '\n';
        return 1;
    }
}
