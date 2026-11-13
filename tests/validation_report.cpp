#include "reference/UniformFieldSolution.hpp"
#include "aurora/physics/Gyromotion.hpp"
#include "aurora/physics/Constants.hpp"
#include "aurora/physics/Relativity.hpp"
#include "aurora/simulation/TrajectoryRunner.hpp"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <sstream>
#include <string>

namespace {
constexpr double tau = 2.0 * std::numbers::pi;

struct Measurement {
    double dt{}, position_error{}, momentum_error{}, energy_error{}, phase_error{}, period_error{},
        radius_error{};
};

Measurement measure(std::size_t per_period, std::size_t periods) {
    using namespace aurora;
    const Vec3 b{0, 0, 3e-5};
    const UniformField field(b);
    const ParticleState initial{
        {0, 0, 0}, momentum_from_energy(constants::eV_to_J(10000), {std::sqrt(3.0), 0, 1}), 0};
    const double period = gyroperiod_s(initial.momentum_kg_m_s, b);
    const double radius = larmor_radius_m(initial.momentum_kg_m_s, b);
    Measurement m;
    m.dt = period / static_cast<double>(per_period);
    double previous_phase = 0, accumulated_phase = 0;
    Vec3 quarter{}, half{};
    const auto summary = run_trajectory(initial, field, {m.dt, per_period * periods, 1},
        [&](std::size_t step, const ParticleState& s) {
            const auto exact = test_reference::uniform_field_solution(initial, b, s.time_s);
            m.position_error = std::max(m.position_error, norm(s.position_m - exact.position_m));
            m.momentum_error = std::max(m.momentum_error,
                norm(s.momentum_kg_m_s - exact.momentum_kg_m_s) / norm(initial.momentum_kg_m_s));
            const double phase = std::atan2(s.momentum_kg_m_s.y, s.momentum_kg_m_s.x);
            accumulated_phase += std::remainder(phase - previous_phase, tau);
            previous_phase = phase;
            if (step == per_period / 4)
                quarter = s.position_m;
            if (step == per_period / 2)
                half = s.position_m;
        });
    m.energy_error = summary.max_relative_energy_error;
    m.phase_error = std::abs(accumulated_phase - tau * static_cast<double>(periods));
    const double measured_period = tau * summary.final_state.time_s / accumulated_phase;
    m.period_error = std::abs(measured_period / period - 1.0);
    // Circumradius of three projected trajectory points, independent of the pusher formula.
    quarter.z = 0;
    half.z = 0;
    const double area_twice = norm(cross(quarter, half));
    const double measured_radius =
        norm(quarter) * norm(half) * norm(half - quarter) / (2.0 * area_twice);
    m.radius_error = std::abs(measured_radius / radius - 1.0);
    return m;
}
}

int main(int argc, char** argv) {
    try {
        std::filesystem::path output = "results/week1_validation.md";
        if (argc == 3 && std::string(argv[1]) == "--output")
            output = argv[2];
        else if (argc != 1)
            throw std::invalid_argument("Usage: aurora_validation_report [--output report.md]");
        std::ostringstream report;
        report
            << std::setprecision(10) << "# Проверка первой недели Aurora\n\n"
            << "Отчёт автоматически сформирован `aurora_validation_report` (Aurora 0.1.0).\n\n"
            << "Компилятор: " <<
#ifdef __VERSION__
            __VERSION__
#else
            "не указан"
#endif
            << ". Режим: " <<
#ifdef NDEBUG
            "NDEBUG (оптимизированная сборка)"
#else
            "с проверками assert"
#endif
            << ".\n\n"
            << "Условия: электрон, K = 10 кэВ, B = (0, 0, 30 мкТл), питч-угол 60°, E = 0. "
            << "Используется double, фиксированный шаг и синхронизация импульса для диагностики.\n\n"
            << "## Сходимость за один аналитический гиропериод\n\n"
            << "Ошибки положения и импульса — максимумы по всем шагам; ошибка импульса нормирована на |p₀|. "
            << "Период измерен по развёрнутой фазе импульса; радиус — по трём точкам проекции траектории.\n\n"
            << "| Шагов/период | dt, с | Ошибка r, м | Ошибка p / p₀ | Ошибка фазы, рад | Отн. ошибка T | Отн. ошибка r_L | Отн. ошибка энергии | Порядок r | Порядок p |\n"
            << "|---|---|---|---|---|---|---|---|---|---|\n";
        bool passed = true;
        Measurement previous{};
        for (const std::size_t n : {40u, 80u, 160u, 320u}) {
            const auto m = measure(n, 1);
            const double order = previous.position_error > 0
                ? std::log2(previous.position_error / m.position_error)
                : 0;
            report << '|' << n << '|' << m.dt << '|' << m.position_error << '|' << m.momentum_error
                   << '|' << m.phase_error << '|' << m.period_error << '|' << m.radius_error << '|'
                   << m.energy_error << '|';
            if (previous.position_error > 0) {
                const double momentum_order = std::log2(previous.momentum_error / m.momentum_error);
                report << order << '|' << momentum_order;
                passed = passed && order > 1.9 && order < 2.1;
                passed = passed && momentum_order > 1.9 && momentum_order < 2.1;
            } else
                report << "—|—";
            report << "|\n";
            passed = passed && std::isfinite(m.radius_error) && m.energy_error < 1e-9;
            previous = m;
        }
        const auto long_run = measure(200, 100);
        passed = passed && long_run.energy_error < 1e-9 && long_run.phase_error > 0 &&
            long_run.phase_error < 0.06 && long_run.period_error < 1e-4 &&
            long_run.radius_error < 2e-4;
        report
            << "\n## Длительный расчёт\n\n"
            << "100 гиропериодов, 200 шагов на период, всего 20 000 шагов.\n\n"
            << "- Максимальное относительное отклонение энергии: " << long_run.energy_error << ".\n"
            << "- Ошибка фазы в конце: " << long_run.phase_error << " рад.\n"
            << "- Относительная ошибка периода: " << long_run.period_error << ".\n"
            << "- Относительная ошибка радиуса: " << long_run.radius_error << ".\n"
            << "- Максимальная ошибка положения: " << long_run.position_error << " м.\n\n"
            << "## Критерии\n\n"
            << "Порядок ошибок положения и импульса: 1.9–2.1. Отклонение энергии: < 10⁻⁹. "
            << "В длительном расчёте ошибка периода < 10⁻⁴, радиуса < 2·10⁻⁴, фазы < 0.06 рад.\n\n"
            << "Статус: **" << (passed ? "PASS" : "FAIL") << "**.\n\n"
            << "Сохранение энергии не устраняет фазовую ошибку. Все сравнения выполняются с непрерывным "
            << "аналитическим решением; сам интегратор его не использует. Нулевое поле, произвольная "
            << "ориентация B, релятивистские преобразования и ошибочный ввод проверяются отдельно через CTest.\n";
        if (output.has_parent_path())
            std::filesystem::create_directories(output.parent_path());
        std::ofstream file(output);
        file.exceptions(std::ios::failbit | std::ios::badbit);
        file << report.str();
        file.close();
        std::cout << (passed ? "PASS: " : "FAIL: ") << output.string() << '\n';
        return passed ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "aurora_validation_report: " << e.what() << '\n';
        return 1;
    }
}
