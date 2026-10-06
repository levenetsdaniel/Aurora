#!/usr/bin/env python3
"""Build one illustrated Markdown report per implemented physics experiment."""

import atexit
import csv
import math
import os
from pathlib import Path
import sys
import tempfile


def rows(directory, filename):
    with (directory / filename).open(newline="", encoding="utf-8") as stream:
        result = list(csv.DictReader(stream))
    if not result:
        raise ValueError(f"Empty experiment data: {filename}")
    return result


def value(row, key):
    return float(row[key])


def write_report(directory, filename, content):
    (directory / filename).write_text(content.strip() + "\n", encoding="utf-8")


def save_plot(figure, directory, filename, plt):
    figure.savefig(directory / filename, dpi=170, bbox_inches="tight")
    plt.close(figure)


def main():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python3 src/experiments/make_reports.py RESULT_DIRECTORY")
    directory = Path(sys.argv[1]).resolve()
    reports = directory / "reports"
    reports.mkdir(parents=True, exist_ok=True)
    cache = tempfile.TemporaryDirectory(prefix="aurora-matplotlib-")
    atexit.register(cache.cleanup)
    os.environ["MPLCONFIGDIR"] = cache.name
    os.environ["XDG_CACHE_HOME"] = cache.name
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.grid": True,
        "grid.alpha": 0.28,
        "figure.facecolor": "white",
    })

    orbit = sorted(rows(directory, "uniform_orbit.csv"),
                   key=lambda row: value(row, "dt_s"))
    event = sorted(rows(directory, "atmosphere_event.csv"),
                   key=lambda row: value(row, "dt_s"))
    energy = rows(directory, "energy_trace.csv")
    mirror = rows(directory, "mirror_trace.csv")
    mirror_summary = rows(directory, "mirror_summary.csv")[0]
    loss_cone = rows(directory, "loss_cone.csv")

    # 01: Relativistic Boris motion in a uniform field and energy conservation.
    dt = [value(row, "dt_s") for row in orbit]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    axes[0].loglog(dt, [value(row, "position_error_m") for row in orbit],
                   "o-", color="#2467a2")
    axes[0].set(xlabel="Шаг dt, с", ylabel="Ошибка координаты, м",
                title="Траектория после четырёх оборотов")
    axes[1].loglog(dt, [abs(value(row, "measured_period_s") /
                            value(row, "exact_period_s") - 1) for row in orbit],
                   "o-", color="#a35122")
    axes[1].set(xlabel="Шаг dt, с", ylabel="Относительная ошибка периода",
                title="Сходимость гиропериода")
    fig.tight_layout()
    save_plot(fig, reports, "01_orbit_convergence.png", plt)

    fig, ax = plt.subplots(figsize=(8, 4))
    for field, color in (("uniform", "#2467a2"), ("dipole", "#a35122")):
        field_rows = [row for row in energy if row["field"] == field]
        ax.plot([value(row, "time_s") * 1e6 for row in field_rows],
                [value(row, "relative_energy_error") for row in field_rows],
                label=field, color=color)
    ax.set(xlabel="Время, мкс", ylabel="(E − E₀) / E₀",
           title="Энергия в статическом магнитном поле")
    ax.legend()
    fig.tight_layout()
    save_plot(fig, reports, "01_energy.png", plt)

    coarse_orbit, fine_orbit = orbit[-1], orbit[0]
    max_energy_error = max(abs(value(row, "relative_energy_error")) for row in energy)
    orbit_table = "\n".join(
        f"| {row['steps_per_period']} | {value(row, 'position_error_m'):.4g} | "
        f"{abs(value(row, 'measured_period_s') / value(row, 'exact_period_s') - 1):.4g} |"
        for row in sorted(orbit, key=lambda row: int(row["steps_per_period"])))
    write_report(reports, "01_integrator.md", f"""
# Эксперимент 1. Движение электрона и проверка интегратора

## Физический вопрос

Воспроизводит ли интегратор Бориса винтовое движение электрона в однородном
магнитном поле, и сохраняется ли энергия, когда электрического поля нет?

## Постановка

Электрон 10 кэВ движется в поле 50 мкТл под наклоном к `B`. Аналитические
радиус и период: `rL = p⊥/(|q|B)` и `T = 2πγm/(|q|B)`. Расчёт продолжается
четыре периода с шагами `T/50`, `T/100`, `T/200`, `T/400`. Отдельные трассы
проверяют энергию в однородном поле и земном диполе.

![Сходимость координаты и периода](01_orbit_convergence.png)

| Шагов за период | Ошибка координаты, м | Относительная ошибка периода |
| ---: | ---: | ---: |
{orbit_table}

![Сохранение энергии](01_energy.png)

## Вывод и ограничения

При уменьшении шага ошибка координаты после четырёх оборотов снизилась с
**{value(coarse_orbit, 'position_error_m'):.3g} до {value(fine_orbit, 'position_error_m'):.3g} м**.
Максимальная относительная ошибка
энергии в сохранённых трассах — **{max_energy_error:.3g}**. Это численная
проверка выбранных случаев, а не гарантия точности любой траектории.
Одна лишь сохранность энергии не доказывает правильность фазы или координаты.

Исходные данные: [орбита](../uniform_orbit.csv),
[энергия](../energy_trace.csv). Код: `src/experiments/TrajectoryExperiments.hpp`.
[Как воспроизвести](../../../src/experiments/README.md).
""")

    # 02: First crossing of a spherical atmospheric hand-off boundary.
    event_dt = [value(row, "dt_s") for row in event]
    fig, ax = plt.subplots(figsize=(8, 4))
    time_line, = ax.loglog(event_dt,
        [value(row, "time_error_s") for row in event],
        "o-", color="#326f51", label="Ошибка времени, с")
    residual_axis = ax.twinx()
    residual_line, = residual_axis.loglog(event_dt,
        [value(row, "sphere_residual_m") for row in event],
        "s--", color="#8b4c91", label="Остаток сферы, м")
    ax.set(xlabel="Шаг dt, с", ylabel="Ошибка времени, с",
           title="Первое пересечение условной сферы 150 км")
    residual_axis.set_ylabel("Остаток сферы, м")
    ax.legend([time_line, residual_line],
              [time_line.get_label(), residual_line.get_label()], loc="upper left")
    fig.tight_layout()
    save_plot(fig, reports, "02_event_location.png", plt)

    coarse_event, fine_event = event[-1], event[0]
    event_table = "\n".join(
        f"| {row['steps_per_period']} | {value(row, 'time_error_s'):.4g} | "
        f"{value(row, 'sphere_residual_m'):.4g} |"
        for row in sorted(event, key=lambda row: int(row["steps_per_period"])))
    write_report(reports, "02_event_location.md", f"""
# Эксперимент 2. Локализация входа в атмосферу

## Физический вопрос

Насколько точно движок определяет **первое** пересечение атмосферной границы,
а не просто последний шаг за ней?

## Постановка

Электрон 10 кэВ в однородном поле 50 мкТл проходит винтовую траекторию около
сферы радиуса `RЗемли + 150 км`. Её аналитическая траектория позволяет найти
эталонное время пересечения. Сравниваются четыре шага интегрирования.
Граница 150 км — **условие передачи** электрона атмосферной модели, а не
место столкновения с Землёй.

![Ошибки времени и положения события](02_event_location.png)

| Шагов за период | Ошибка времени, с | Остаток на сфере, м |
| ---: | ---: | ---: |
{event_table}

## Вывод и ограничения

Ошибка времени снизилась с **{value(coarse_event, 'time_error_s'):.3g} до {value(fine_event, 'time_error_s'):.3g} с**
при уменьшении шага. Остаток
`||rсобытия| − Rграницы|` измеряется **в метрах** и не равен ошибке времени
или полной ошибке траектории. Опыт проверяет локатор на аналитическом
однородном поле; в неоднородном поле нужна отдельная проверка сходимости.

Исходные данные: [события](../atmosphere_event.csv).
Код: `src/simulation/EventDetector.hpp` и `src/experiments/TrajectoryExperiments.hpp`.
[Как воспроизвести](../../../src/experiments/README.md).
""")

    # 03: One dipole trajectory reaching its first magnetic mirror.
    first_b = value(mirror[0], "B_T")
    first_mu = value(mirror[0], "mu_J_T")
    first_p = abs(value(mirror[0], "p_parallel_kg_m_s"))
    time_ms = [value(row, "time_s") * 1e3 for row in mirror]
    fig, axes = plt.subplots(3, 1, figsize=(9, 7), sharex=True)
    axes[0].plot(time_ms, [value(row, "B_T") / first_b for row in mirror],
                 color="#2467a2")
    axes[0].set(ylabel="B / B₀", title="Электрон приближается к зеркалу")
    axes[1].plot(time_ms, [value(row, "p_parallel_kg_m_s") / first_p
                           for row in mirror], color="#7b3e92")
    axes[1].axhline(0, color="black", linewidth=0.7)
    axes[1].set(ylabel="p∥ / |p∥₀|")
    mu_ppm = [(value(row, "mu_J_T") / first_mu - 1.0) * 1e6 for row in mirror]
    axes[2].plot(time_ms, mu_ppm, color="#a35122", alpha=0.55, linewidth=0.55)
    axes[2].set(xlabel="Время, мс", ylabel="(μ / μ₀ − 1), ppm")
    fig.tight_layout()
    save_plot(fig, reports, "03_magnetic_mirror.png", plt)

    mirror_ratio = (value(mirror_summary, "B_event_T") /
                    value(mirror_summary, "B_m_predicted_T"))
    write_report(reports, "03_magnetic_mirror.md", f"""
# Эксперимент 3. Магнитное зеркало

## Физический вопрос

Отражается ли электрон там, где это предсказывает сохранение магнитного
момента `μ = p⊥²/(2mB)`?

## Постановка

Электрон 10 кэВ стартует на высоте 450 км с локальным углом 80° относительно
направления движения вдоль поля к Земле. Поле — статический земной диполь.
Сравниваются поле первого отражения `Bm` и адиабатическая оценка
`Bm,теор = B₀/sin²(80°)`.

![Поле, продольный импульс и магнитный момент](03_magnetic_mirror.png)

## Вывод и ограничения

Исход частицы — **{mirror_summary['outcome']}**. При движении `B` возрастает,
`p∥` достигает нуля; отношение **`Bm/Bm,теор = {mirror_ratio:.7f}`**.
Изменение мгновенного `μ` на графике показано в миллионных долях (`ppm`);
колебания на гиропериоде не следует интерпретировать как вековой дрейф.
Движок останавливается при **первом** зеркале: обратный полёт и период
многократных отражений здесь не моделируются.

Исходные данные: [сводка](../mirror_summary.csv),
[трасса](../mirror_trace.csv). Код: `src/experiments/TrajectoryExperiments.hpp`.
[Как воспроизвести](../../../src/experiments/README.md).
""")

    # 04: Local loss-cone boundary versus dipole L and boundary altitude.
    main_lines = sorted((row for row in loss_cone
                         if value(row, "atmosphere_height_m") == 150000.0 and
                         int(row["steps_per_gyroperiod"]) == 50),
                        key=lambda row: value(row, "L"))
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.errorbar([value(row, "L") for row in main_lines],
                [value(row, "numerical_pitch_deg") for row in main_lines],
                yerr=[(value(row, "mirroring_pitch_deg") -
                       value(row, "precipitating_pitch_deg")) / 2
                      for row in main_lines],
                fmt="o-", capsize=3, color="#2467a2", label="Численный интервал")
    ax.plot([value(row, "L") for row in main_lines],
            [value(row, "theoretical_pitch_deg") for row in main_lines],
            "s--", color="#a35122", label="Адиабатическая оценка")
    ax.set(xlabel="Параметр силовой линии L", ylabel="Локальный угол α₀, градусы",
           title="Граница конуса потерь на высоте 450 км")
    ax.legend()
    fig.tight_layout()
    save_plot(fig, reports, "04_loss_cone_vs_L.png", plt)

    central_colatitude = 0.32 * 180 / math.pi
    central = [row for row in loss_cone
               if abs(value(row, "source_colatitude_deg") - central_colatitude) < 1e-6]
    heights = sorted((row for row in central
                      if int(row["steps_per_gyroperiod"]) == 50),
                     key=lambda row: value(row, "atmosphere_height_m"))
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot([value(row, "atmosphere_height_m") / 1000 for row in heights],
            [value(row, "numerical_pitch_deg") for row in heights],
            "o-", color="#2467a2", label="Расчёт")
    ax.plot([value(row, "atmosphere_height_m") / 1000 for row in heights],
            [value(row, "theoretical_pitch_deg") for row in heights],
            "s--", color="#a35122", label="Теория")
    ax.set(xlabel="Высота условной границы, км", ylabel="Локальный угол α₀, градусы",
           title="Зависимость от определения «попал в атмосферу»")
    ax.legend()
    fig.tight_layout()
    save_plot(fig, reports, "04_boundary_height.png", plt)

    central_150 = [row for row in central
                   if value(row, "atmosphere_height_m") == 150000.0]
    coarse = next(row for row in central_150
                  if int(row["steps_per_gyroperiod"]) == 50)
    fine = next(row for row in central_150
                if int(row["steps_per_gyroperiod"]) == 100)
    max_difference = max(abs(value(row, "numerical_pitch_deg") -
                             value(row, "theoretical_pitch_deg"))
                         for row in main_lines)
    loss_table = "\n".join(
        f"| {value(row, 'L'):.2f} | {value(row, 'numerical_pitch_deg'):.3f} | "
        f"{value(row, 'theoretical_pitch_deg'):.3f} | "
        f"{value(row, 'precipitating_pitch_deg'):.3f}–"
        f"{value(row, 'mirroring_pitch_deg'):.3f} |"
        for row in main_lines)
    write_report(reports, "04_loss_cone.md", f"""
# Эксперимент 4. Граница конуса потерь

## Физический вопрос

Какой локальный угол запуска отделяет высыпающиеся электроны от отражённых,
и зависит ли он от силовой линии и выбранной высоты атмосферной границы?

## Постановка

Электрон 10 кэВ стартует с высоты 450 км в статическом диполе. Гирофаза
одинакова для всех прогонов. Бинарный поиск удерживает одну сторону скобки
`Precipitated`, другую — `Mirrored`; `Escaped` и `Unresolved` вызывают ошибку.
Теоретическая оценка в адиабатическом приближении:
`sin²αc ≈ B(450 км)/B(граница)` вдоль **той же** силовой линии.

![Граница для трёх линий](04_loss_cone_vs_L.png)

| L | αc численно, ° | αc теория, ° | Скобка, ° |
| ---: | ---: | ---: | ---: |
{loss_table}

Усики на графике — половина ширины скобки бинарного поиска, **не**
статистический доверительный интервал.

![Влияние высоты условной границы](04_boundary_height.png)

## Вывод и ограничения

На трёх линиях максимальное расхождение численной и теоретической середины
скобки — **{max_difference:.3f}°**; ширина скобки менее 0,05°.
При переходе от `T/50` к `T/100` центральная оценка изменилась на
{abs(value(coarse, 'numerical_pitch_deg') - value(fine, 'numerical_pitch_deg')):.4f}° —
неразличимо при текущем разрешении поиска. Это **локальный угол на
высоте 450 км**, не экваториальный угол. Перенос условной границы со 100 на
200 км меняет критический угол: определение «высыпался» существенно.

Здесь нет статистического источника, возмущённого поля или реальных
спутниковых данных. Согласие подтверждает модельный механизм магнитного
зеркала, но не является прогнозом сияния.

Исходные данные: [полный CSV](../loss_cone.csv).
Код: `src/experiments/LossConeExperiments.hpp`.
[Как воспроизвести](../../../src/experiments/README.md).
""")

    write_report(reports, "README.md", """
# Отчёты по физическим экспериментам

Каждый реализованный опыт имеет отдельный отчёт с постановкой, графиками,
численными результатами, физическим выводом и ограничениями:

1. [Движение электрона и проверка интегратора](01_integrator.md)
2. [Локализация входа в атмосферу](02_event_location.md)
3. [Магнитное зеркало](03_magnetic_mirror.md)
4. [Граница конуса потерь](04_loss_cone.md)

Все рисунки лежат рядом с отчётами. Исходные CSV — в родительском каталоге.
Графики и отчёты можно [пересоздать по инструкции](../../../src/experiments/README.md).
""")
    print(f"Wrote four Markdown reports and six figures to {reports}")


if __name__ == "__main__":
    main()
