#include "aurora/fields/DipoleField.hpp"
#include "aurora/physics/ParticleProperties.hpp"
#include "aurora/io/TrajectoryCsv.hpp"
#include "aurora/physics/Relativity.hpp"
#include "aurora/simulation/TrajectoryRunner.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>

int main() {
    try {
        // Учебный пример: электрон стартует на магнитном экваторе и отражается в диполе.
        const double energy_keV = 100.0;
        const double pitch_deg = 45.0;
        const double equatorial_distance_RE = 2.0;
        const double dt_s = 2e-8;
        const double max_time_s = 0.1;

        const aurora::DipoleField field; // Момент вдоль -z; на экваторе B направлено вдоль +z.
        const double pitch = pitch_deg * std::numbers::pi / 180.0;
        if (pitch_deg<=0 || pitch_deg>=90 || equatorial_distance_RE<=1.02)
            throw std::invalid_argument("Для примера задайте 0<pitch_deg<90 и расстояние>1.02 R_E");
        const aurora::ParticleState electron{
            {equatorial_distance_RE * aurora::constants::earth_radius_m, 0, 0},
            aurora::momentum_from_energy(aurora::constants::eV_to_J(energy_keV*1000),
                {std::sin(pitch), 0, std::cos(pitch)}), 0
        };
        aurora::RunConfig config{dt_s, 10000000, 100};
        config.max_duration_s = max_time_s;
        config.inner_radius_m = 1.02 * aurora::constants::earth_radius_m;
        config.outer_radius_m = 4.0 * aurora::constants::earth_radius_m;
        config.detect_mirrors = true;
        config.stop_at_first_mirror = true;

        const auto output = std::filesystem::absolute("results");
        std::filesystem::create_directories(output);
        std::ofstream csv(output / "trajectory.csv");
        csv.exceptions(std::ios::failbit | std::ios::badbit);
        aurora::write_magnetic_trajectory_header(csv);
        const auto started = std::chrono::steady_clock::now();
        const auto result = aurora::run_trajectory(electron, field, config,
            [&](std::size_t step, const aurora::ParticleState& state) {
                aurora::write_magnetic_trajectory_row(csv, step, state, field);
            });
        csv.close();
        const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();

        std::ofstream events(output / "events.csv");
        events.exceptions(std::ios::failbit | std::ios::badbit);
        events << "type,time_s,x_m,y_m,z_m,B_T,bracket_width_s\n" << std::setprecision(17);
        for (const auto& e : result.events) {
            const char* name = e.type==aurora::EventType::Mirror ? "mirror"
                : (e.type==aurora::EventType::InnerBoundary ? "inner_boundary" : "outer_boundary");
            events << name << ',' << e.state.time_s << ',' << e.state.position_m.x << ','
                   << e.state.position_m.y << ',' << e.state.position_m.z << ','
                   << aurora::norm(field.at(e.state.position_m,e.state.time_s)) << ',' << e.bracket_width_s << '\n';
        }
        events.close();
        std::cout << std::setprecision(7)
                  << "Aurora: магнитное отражение электрона в дипольном поле\n"
                  << "Энергия: " << energy_keV << " кэВ; начальный питч-угол: " << pitch_deg << " градусов\n"
                  << "Причина остановки: " << aurora::stop_reason_name(result.stop_reason) << '\n'
                  << "Физическое время: " << result.final_state.time_s << " с; шагов: " << result.completed_steps << '\n'
                  << "Отклонение энергии: " << result.max_relative_energy_error << '\n'
                  << "Максимальное отклонение первого инварианта mu: " << result.max_relative_magnetic_moment_error << '\n';
        if (result.stop_reason==aurora::StopReason::Mirror) {
            const double expected = aurora::predicted_mirror_field_T(electron,field);
            const double measured = aurora::norm(field.at(result.final_state.position_m,result.final_state.time_s));
            std::cout << "B отражения: " << measured << " Тл; адиабатическая оценка: " << expected << " Тл\n"
                      << "Относительное отличие: " << std::abs(measured/expected-1) << '\n';
        }
        std::cout << "Время вычислений: " << elapsed << " с\n"
                  << "Траектория: " << (output / "trajectory.csv").string() << '\n'
                  << "События: " << (output / "events.csv").string() << '\n';
        if (result.stop_reason==aurora::StopReason::NumericalError)
            throw std::runtime_error(result.error_message);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Ошибка: " << error.what() << '\n';
        return 1;
    }
}
