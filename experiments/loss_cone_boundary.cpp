#include "aurora/io/ExperimentFiles.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <vector>

namespace {
using namespace aurora;

constexpr double angular_tolerance_deg = 0.0001;
constexpr std::array scan_angles_deg{5., 10., 15., 20., 30., 45., 60., 75., 85.};

struct Boundary {
    int direction;
    double dt_s;
    double entry_angle_deg;
    double reflection_angle_deg;
};

std::filesystem::path new_session(const std::filesystem::path& root) {
    std::filesystem::create_directories(root);
    const auto stamp = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::system_clock::now().time_since_epoch())
                           .count();

    for (int index = 0; index < 10000; ++index) {
        const auto path = std::filesystem::absolute(
            root / ("loss-cone-" + std::to_string(stamp) + "-" + std::to_string(index)));
        if (std::filesystem::create_directory(path))
            return path;
    }
    throw std::runtime_error("Cannot create loss cone experiment directory");
}

class Experiment {
  public:
    explicit Experiment(const std::filesystem::path& directory) : directory_(directory) {
        base_.max_time_s = 0.5;
        base_.max_steps = 100000000;
        base_.sample_every = 10000;
        base_.stop_at_first_mirror = true;

        std::ofstream parameters(directory_ / "base_parameters.txt");
        parameters.exceptions(std::ios::failbit | std::ios::badbit);
        write_experiment_parameters(parameters, base_);
        parameters.close();

        csv_.exceptions(std::ios::failbit | std::ios::badbit);
        csv_.open(directory_ / "cases.csv");
        csv_ << std::setprecision(17)
             << "case_id,stage,direction,dt_s,local_angle_deg,pitch_deg,energy_keV,gyro_phase_deg,"
                "outcome,time_s,altitude_m,max_energy_error,max_mu_error,run_directory\n";
    }

    Boundary search(int direction, double dt_s) {
        double entry = 0, reflection = 0;
        bool found_entry = false, found_reflection = false;

        for (double angle : scan_angles_deg) {
            const auto outcome = observe(angle, direction, dt_s, "scan");
            if (outcome == ExperimentOutcome::AtmosphereReached) {
                if (found_reflection)
                    throw std::runtime_error("Coarse scan contains more than one transition");
                entry = angle;
                found_entry = true;
            } else if (!found_reflection) {
                reflection = angle;
                found_reflection = true;
            }
        }
        if (!found_entry || !found_reflection || entry >= reflection)
            throw std::runtime_error("Coarse scan did not bracket entry and reflection");

        while (reflection - entry > angular_tolerance_deg) {
            const double middle = entry + (reflection - entry) / 2;
            if (observe(middle, direction, dt_s, "bisection") ==
                ExperimentOutcome::AtmosphereReached)
                entry = middle;
            else
                reflection = middle;
        }
        return {direction, dt_s, entry, reflection};
    }

    double analytical_angle_deg() const {
        return prepare_experiment(base_).loss_cone.critical_angle_rad * 180 / std::numbers::pi;
    }

  private:
    ExperimentOutcome observe(double local_angle, int direction, double dt_s, const char* stage) {
        auto parameters = base_;
        parameters.dt_s = dt_s;
        parameters.pitch_deg = direction > 0 ? local_angle : 180 - local_angle;
        const auto saved = save_experiment(parameters, directory_ / "cases");
        const auto& run = saved.result.trajectory;
        const auto outcome = saved.result.outcome;
        const double altitude = norm(run.final_state.position_m) - constants::earth_radius_m;

        csv_ << case_id_++ << ',' << stage << ',' << direction << ',' << dt_s << ',' << local_angle
             << ',' << parameters.pitch_deg << ',' << parameters.energy_keV << ','
             << parameters.gyro_phase_deg << ',' << outcome_name(outcome) << ','
             << run.final_state.time_s << ',' << altitude << ',' << run.max_relative_energy_error
             << ',' << run.max_relative_magnetic_moment_error << ','
             << std::filesystem::relative(saved.directory, directory_).generic_string() << '\n';
        csv_.flush();
        std::cout << "Case " << case_id_ << ": " << (direction > 0 ? "+B" : "-B") << ", dt=" << dt_s
                  << ", angle=" << std::setprecision(10) << local_angle << " -> "
                  << outcome_name(outcome) << '\n'
                  << std::flush;

        if (outcome != ExperimentOutcome::AtmosphereReached &&
            outcome != ExperimentOutcome::Reflected)
            throw std::runtime_error("Boundary search requires a resolved entry or reflection: " +
                saved.directory.string());
        if (!std::isfinite(run.max_relative_energy_error) || run.max_relative_energy_error >= 1e-9)
            throw std::runtime_error(
                "Energy conservation check failed: " + saved.directory.string());
        return outcome;
    }

    ExperimentParameters base_;
    std::filesystem::path directory_;
    std::ofstream csv_;
    std::size_t case_id_{};
};
}

int main(int argc, char** argv) {
    try {
        std::filesystem::path root = "results/loss-cone";
        if (argc == 3 && std::string(argv[1]) == "--output")
            root = argv[2];
        else if (argc != 1)
            throw std::invalid_argument("Usage: aurora_loss_cone_boundary [--output directory]");

        const auto directory = new_session(root);
        std::cout << "Results: " << directory.string() << '\n' << std::flush;
        Experiment experiment(directory);
        std::vector<Boundary> boundaries;
        for (double dt : {4e-8, 2e-8, 1e-8})
            for (int direction : {1, -1})
                boundaries.push_back(experiment.search(direction, dt));

        const double analytical = experiment.analytical_angle_deg();
        std::ofstream comparison(directory / "comparison.csv");
        comparison.exceptions(std::ios::failbit | std::ios::badbit);
        comparison << std::setprecision(17)
                   << "direction,dt_s,entry_angle_deg,reflection_angle_deg,midpoint_deg,"
                      "half_width_deg,analytical_deg,offset_deg,relative_offset_percent\n";
        for (const auto& boundary : boundaries) {
            const double midpoint = (boundary.entry_angle_deg + boundary.reflection_angle_deg) / 2;
            const double half_width =
                (boundary.reflection_angle_deg - boundary.entry_angle_deg) / 2;
            comparison << boundary.direction << ',' << boundary.dt_s << ','
                       << boundary.entry_angle_deg << ',' << boundary.reflection_angle_deg << ','
                       << midpoint << ',' << half_width << ',' << analytical << ','
                       << midpoint - analytical << ',' << 100 * (midpoint - analytical) / analytical
                       << '\n';
            std::cout << std::setprecision(12) << (boundary.direction > 0 ? "+B" : "-B")
                      << ", dt=" << boundary.dt_s << ": [" << boundary.entry_angle_deg << ", "
                      << boundary.reflection_angle_deg << "] deg; analytical=" << analytical
                      << '\n';
        }
        comparison.close();

        std::cout << "PASS: " << directory.string() << '\n';
        const double finest_step = std::min_element(
            boundaries.begin(), boundaries.end(), [](const Boundary& a, const Boundary& b) {
                return a.dt_s < b.dt_s;
            })->dt_s;
        std::cout << "Итоговый интервал границы при шаге " << std::scientific
                  << std::setprecision(0) << finest_step << " с (угол к направлению поля):\n";
        for (const auto& boundary : boundaries) {
            if (boundary.dt_s == finest_step)
                std::cout << std::fixed << std::setprecision(8)
                          << (boundary.direction > 0 ? "  +B: [" : "  -B: [")
                          << boundary.entry_angle_deg << "; " << boundary.reflection_angle_deg
                          << "] градусов\n";
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Loss cone experiment: " << error.what() << '\n';
        return 1;
    }
}
