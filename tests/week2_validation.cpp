#include "aurora/fields/DipoleField.hpp"
#include "aurora/physics/ParticleProperties.hpp"
#include "aurora/physics/Relativity.hpp"
#include "aurora/simulation/TrajectoryRunner.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <sstream>
#include <vector>

namespace {
using namespace aurora;
struct Measurement {
    double energy_keV, pitch_deg, dt_s, seconds;
    ParticleState initial;
    RunSummary run;
    double mirror_error, latitude_deg, altitude_km;
};
Measurement measure(double keV,double pitch_deg,double dt) {
    const DipoleField field;
    const double pitch=pitch_deg*std::numbers::pi/180;
    const ParticleState initial{{2*constants::earth_radius_m,0,0},
        momentum_from_energy(constants::eV_to_J(keV*1000),{std::sin(pitch),0,std::cos(pitch)}),0};
    RunConfig c{dt,100000000,1000};
    c.max_duration_s=1; c.detect_mirrors=true; c.stop_at_first_mirror=true;
    c.inner_radius_m=1.02*constants::earth_radius_m; c.outer_radius_m=4*constants::earth_radius_m;
    c.event_tolerance_s=1e-13;
    const auto start=std::chrono::steady_clock::now();
    auto run=run_trajectory(initial,field,c);
    const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    if (run.stop_reason!=StopReason::Mirror) throw std::runtime_error("No mirror: "+std::string(stop_reason_name(run.stop_reason))+" "+run.error_message);
    const Vec3 p=run.final_state.position_m;
    const double error=std::abs(norm(field.at(p,0))/predicted_mirror_field_T(initial,field)-1);
    return {keV,pitch_deg,dt,elapsed,initial,std::move(run),error,
        std::asin(p.z/norm(p))*180/std::numbers::pi,(norm(p)-constants::earth_radius_m)/1000};
}
}
int main(int argc,char** argv) {
    try {
        std::filesystem::path output="results/week2_validation.md";
        if (argc==3 && std::string(argv[1])=="--output") output=argv[2];
        else if (argc!=1) throw std::invalid_argument("Usage: aurora_week2_report [--output report.md]");
        const DipoleField field;
        bool passed=true;
        std::ostringstream out;
        out << std::setprecision(10)
            << "# Aurora: результаты второй недели\n\n"
            << "Автоматический отчёт `aurora_week2_report`. Компилятор: "
#ifdef __VERSION__
            << __VERSION__
#else
            << "unknown"
#endif
#ifdef NDEBUG
            << "; NDEBUG (Release).\n\n"
#else
            << "; Debug.\n\n"
#endif
            << "Диполь: R_E=6371000 м, B_eq(R_E)=31.2 мкТл, ось момента (0,0,−1). "
            << "Старт (2 R_E,0,0), начальный импульс в плоскости xz, E=0. "
            << "Остановка на первом отражении, допуск локализации 10⁻¹³ с.\n\n"
            << "## Энергия и питч-угол\n\n"
            << "Фиксированный шаг 4·10⁻⁸ с. Положение отражения задано магнитной широтой и высотой над сферической Землёй. "
            << "Отклонение μ — максимум мгновенного |μ/μ₀−1| по всем шагам, включая гироколебания.\n\n"
            << "| K, кэВ | α₀, град | t отражения, с | Широта, град | Высота, км | B отражения, Тл | B₀/sin²α₀, Тл | Отн. отличие B | max Δμ/μ₀ | max ΔK/K₀ |\n"
            << "|---|---|---|---|---|---|---|---|---|---|\n";
        std::vector<Measurement> grid;
        for (double energy : {10.,100.,1000.}) for (double angle : {30.,45.,60.}) {
            auto m=measure(energy,angle,4e-8);
            const double measured=norm(field.at(m.run.final_state.position_m,0));
            out << '|' << energy << '|' << angle << '|' << m.run.final_state.time_s << '|'
                << m.latitude_deg << '|' << m.altitude_km << '|' << measured << '|'
                << predicted_mirror_field_T(m.initial,field) << '|' << m.mirror_error << '|'
                << m.run.max_relative_magnetic_moment_error << '|' << m.run.max_relative_energy_error << "|\n";
            passed=passed && m.mirror_error<0.01 && m.run.max_relative_energy_error<1e-9
                && m.run.max_relative_magnetic_moment_error<0.01;
            grid.push_back(std::move(m));
        }
        out << "\nПри заданном питч-угле поле отражения близко к одинаковому для разных энергий, "
            << "а время достижения отражения уменьшается при увеличении энергии. "
            << "Это ожидаемо для адиабатического движения в стационарном поле; конечный ларморовский радиус даёт отличия.\n\n"
            << "## Сходимость точки и времени отражения\n\n"
            << "На длинной гиротраектории грубый шаг может давать немонотонную ошибку точки отражения "
            << "из-за накопления фазы. Поэтому здесь используется диапазон меньших шагов, а не один "
            << "расчёт с хорошим сохранением энергии.\n\n"
            << "K=1000 кэВ, α₀=45°. Эталон — отдельный полный гирорасчёт с dt=2.5·10⁻⁹ с, "
            << "а не адиабатическая формула. Ошибка положения измеряется по трём координатам в метрах.\n\n"
            << "| dt, с | t отражения, с | Отличие времени от эталона, с | Отличие положения от эталона, м | Ширина скобки, с |\n"
            << "|---|---|---|---|---|\n";
        const auto reference=measure(1000,45,2.5e-9);
        double previous_time=0,previous_position=0;
        for (double dt : {4e-8,2e-8,1e-8,5e-9}) {
            const auto m=measure(1000,45,dt);
            const double time_error=std::abs(m.run.final_state.time_s-reference.run.final_state.time_s);
            const double position_error=norm(m.run.final_state.position_m-reference.run.final_state.position_m);
            out << '|' << dt << '|' << m.run.final_state.time_s << '|' << time_error << '|'
                << position_error << '|' << m.run.events[0].bracket_width_s << "|\n";
            if (previous_time>0) passed=passed && time_error<previous_time && position_error<previous_position;
            previous_time=time_error; previous_position=position_error;
        }
        const auto& benchmark=grid[4];
        const double throughput=static_cast<double>(benchmark.run.completed_steps)/benchmark.seconds;
        out << "\nДополнительно CTest проверяет второй порядок времени пересечения сферы по точному решению "
            << "в однородном поле. Уменьшение скобки локализатора отдельно от dt не устраняет ошибку траектории.\n\n"
            << "## Стоимость полного гирорасчёта\n\n"
            << "Один запуск 100 кэВ, 45°, dt=4·10⁻⁸ с до отражения, с диагностикой каждого шага, без CSV:\n\n"
            << "- Шагов: " << benchmark.run.completed_steps << ".\n"
            << "- Физическое время: " << benchmark.run.final_state.time_s << " с.\n"
            << "- Время вычислений: " << benchmark.seconds << " с.\n"
            << "- Производительность: " << throughput << " шагов/с.\n"
            << "- Линейная оценка для 10 000 таких частиц последовательно: " << 10000*benchmark.seconds/3600 << " ч.\n\n"
            << "Это однократное измерение на текущем компьютере, не универсальный benchmark. "
            << "Оценка ансамбля не включает вывод, разброс траекторий и параллелизм. "
            << "Полный гирорасчёт дорог; до ансамбля необходимо выбирать число частиц и шаг по исследованиям сходимости.\n\n"
            << "## Критерии и ограничения\n\n"
            << "Для таблицы: отличие B отражения от адиабатической оценки <1%, мгновенное отклонение μ <1%, "
            << "отклонение энергии <10⁻⁹. При уменьшении dt ошибки времени и положения отражения должны уменьшаться.\n\n"
            << "Статус: **" << (passed ? "PASS" : "FAIL") << "**.\n\n"
            << "μ=p_perp²/(2mB) — диагностический адиабатический инвариант, не строго сохраняющаяся величина "
            << "полной гиротраектории. Численная сходимость не устраняет физические гироколебания μ. "
            << "Диполь — идеализированное заданное поле; атмосфера, столкновения и конус потерь в этом отчёте не моделируются.\n\n"
            << "Определение отражения и адиабатическая оценка: "
            << "[лекции R. Fitzpatrick, Magnetic Mirrors](https://farside.ph.utexas.edu/teaching/plasma/Plasma/node22.html). "
            << "В коде используется релятивистский импульс; при постоянной энергии отношение полей получается "
            << "из сохранения p_perp²/B и |p|.\n";
        if (output.has_parent_path()) std::filesystem::create_directories(output.parent_path());
        std::ofstream file(output); file.exceptions(std::ios::failbit|std::ios::badbit);
        file << out.str(); file.close();
        std::cout << (passed ? "PASS: " : "FAIL: ") << output.string() << '\n';
        return passed ? 0 : 1;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
