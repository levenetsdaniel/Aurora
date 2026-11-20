#include "aurora/io/ExperimentFiles.hpp"
#include "aurora/physics/ParticleProperties.hpp"
#include "aurora/physics/Relativity.hpp"
#include <iomanip>
#include <iostream>
#include <numbers>
#include <cstdlib>
#ifdef AURORA_HAS_VTK
#include "aurora/visualization/VtkViewer.hpp"
#endif

int main() {
    try {
        aurora::ExperimentParameters parameters;
        parameters.energy_keV = 100.0;
        parameters.pitch_deg = 45.0; // 5: вход вдоль B; 175: вход против B.
        parameters.equatorial_distance_RE = 2.0;
        parameters.dt_s = 2e-8;
        parameters.max_time_s = 0.1;
        parameters.atmosphere_height_m = 100000.0;

#ifdef AURORA_HAS_VTK
        if (!std::getenv("AURORA_NO_GUI")) {
            const auto cases =
                aurora::visualization::create_demo_cases(parameters, "results/visualization");
            aurora::visualization::show_demo(cases);
            return 0;
        }
#endif

        const auto prepared = aurora::prepare_experiment(parameters);
        const auto& cone = prepared.loss_cone;

        std::cout << std::setprecision(10)
                  << "Критический питч-угол: " << cone.critical_angle_rad * 180 / std::numbers::pi
                  << " градусов\n"
                  << "Начальный питч-угол: " << parameters.pitch_deg << " градусов\n"
                  << "Аналитический прогноз: "
                  << (prepared.prediction ? aurora::prediction_name(*prepared.prediction)
                                          : "not_applicable")
                  << '\n'
                  << std::flush;

        const auto saved = aurora::save_experiment(parameters, "results");
        const auto& result = saved.result;
        const auto& run = result.trajectory;
        const auto& state = run.final_state;

        std::cout << "Исход: " << aurora::outcome_name(result.outcome)
                  << "; сравнение с прогнозом: " << aurora::agreement_name(result.agreement) << '\n'
                  << "Причина остановки: " << aurora::stop_reason_name(run.stop_reason) << '\n'
                  << "Физическое время: " << state.time_s << " с; шагов: " << run.completed_steps
                  << '\n'
                  << "Отклонение энергии: " << run.max_relative_energy_error << '\n'
                  << "Максимальное отклонение mu: " << run.max_relative_magnetic_moment_error
                  << '\n';

        if (result.outcome == aurora::ExperimentOutcome::AtmosphereReached) {
            std::cout << "Электрон достиг атмосферной границы.\n"
                      << "Координаты входа: " << state.position_m.x << ", " << state.position_m.y
                      << ", " << state.position_m.z << " м\n"
                      << "Высота: "
                      << aurora::norm(state.position_m) - aurora::constants::earth_radius_m
                      << " м\n"
                      << "Энергия входа: "
                      << aurora::constants::J_to_eV(
                             aurora::kinetic_energy_J(state.momentum_kg_m_s)) /
                    1000
                      << " кэВ\n";

        } else if (result.outcome == aurora::ExperimentOutcome::Reflected) {
            const double expected =
                aurora::predicted_mirror_field_T(prepared.initial, prepared.field);
            const double measured = aurora::norm(prepared.field.at(state.position_m, state.time_s));

            std::cout << "B отражения: " << measured << " Тл; адиабатическая оценка: " << expected
                      << " Тл; отличие: " << std::abs(measured / expected - 1) << '\n';
        }

        std::cout << "Время вычислений: " << saved.wall_time_s << " с\n"
                  << "Результаты запуска: " << saved.directory.string() << '\n';

        if (result.outcome == aurora::ExperimentOutcome::NumericalError)
            throw std::runtime_error(run.error_message);

        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Ошибка: " << error.what() << '\n';
        return 1;
    }
}
