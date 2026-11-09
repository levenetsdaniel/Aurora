#include "aurora/diagnostics/UniformFieldSolution.hpp"
#include "aurora/io/TrajectoryCsv.hpp"
#include "aurora/physics/Constants.hpp"
#include "aurora/physics/Relativity.hpp"
#include "aurora/simulation/TrajectoryRunner.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>

int main() {
    try {
        const double energy_keV = 10.0;
        const double pitch_deg = 60.0;
        const double magnetic_field_T = 3e-5; // Поле направлено вдоль +z.
        const std::size_t periods = 10;
        const std::size_t steps_per_period = 200;

        const aurora::UniformField field({0, 0, magnetic_field_T});
        const double pitch = pitch_deg * std::numbers::pi / 180.0;
        const aurora::ParticleState electron{
            {0, 0, 0},
            aurora::momentum_from_energy(
                aurora::constants::eV_to_J(energy_keV * 1000.0),
                {std::sin(pitch), 0, std::cos(pitch)}),
            0.0
        };
        const double period_s = aurora::gyroperiod_s(electron.momentum_kg_m_s, field.value_T());
        if (steps_per_period == 0) throw std::invalid_argument("Число шагов на период должно быть больше нуля");
        if (periods > std::numeric_limits<std::size_t>::max() / steps_per_period)
            throw std::invalid_argument("Слишком большое число шагов");
        const aurora::RunConfig config{period_s / static_cast<double>(steps_per_period),
                                       periods * steps_per_period, 1};

        const auto output = std::filesystem::absolute("results/trajectory.csv");
        std::filesystem::create_directories(output.parent_path());
        std::ofstream csv(output);
        csv.exceptions(std::ios::failbit | std::ios::badbit);
        aurora::write_trajectory_header(csv);
        const auto result = aurora::run_trajectory(electron, field, config,
            [&](std::size_t step, const aurora::ParticleState& state) {
                aurora::write_trajectory_row(csv, step, state);
            });
        csv.close();

        std::cout << std::setprecision(6)
                  << "Aurora: электрон в однородном магнитном поле\n"
                  << "Энергия: " << energy_keV << " кэВ\n"
                  << "Питч-угол: " << pitch_deg << " градусов\n"
                  << "Магнитное поле: " << magnetic_field_T << " Тл\n"
                  << "Рассчитано периодов: " << periods << " (" << result.completed_steps << " шагов)\n"
                  << "Гиропериод: " << period_s << " с\n"
                  << "Ларморовский радиус: "
                  << aurora::larmor_radius_m(electron.momentum_kg_m_s, field.value_T()) << " м\n"
                  << "Максимальное относительное отклонение энергии: " << result.max_relative_energy_error << '\n'
                  << "Траектория сохранена: " << output.string() << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Ошибка: " << error.what() << '\n';
        return 1;
    }
}
