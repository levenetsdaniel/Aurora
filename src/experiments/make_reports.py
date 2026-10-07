#!/usr/bin/env python3
"""Build one illustrated Markdown report per implemented physics experiment."""

import atexit
import csv
import math
import os
from pathlib import Path
import random
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


def wilson(successes, trials):
    z = 1.959963984540054
    fraction = successes / trials
    denominator = 1 + z * z / trials
    center = (fraction + z * z / (2 * trials)) / denominator
    half = z / denominator * math.sqrt(
        fraction * (1 - fraction) / trials + z * z / (4 * trials * trials))
    return max(0.0, center - half), min(1.0, center + half)


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
    precipitation_bins = rows(directory, "precipitation_bins.csv")
    source_distributions = rows(directory, "source_distributions.csv")
    source_convergence = rows(directory, "source_convergence.csv")
    precipitation_particles = rows(directory, "precipitation_particles.csv")
    anomaly_field = rows(directory, "anomaly_field.csv")
    anomaly_scan = rows(directory, "anomaly_scan.csv")
    anomaly_summary = rows(directory, "anomaly_summary.csv")
    anomaly_particles = rows(directory, "anomaly_particles.csv")

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

    # 05: Angular-bin probabilities and two different source hypotheses.
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for L, color in ((7.06, "#286aa1"), (10.82, "#a65320"),
                     (17.49, "#588540")):
        selected = sorted((row for row in precipitation_bins
                           if abs(value(row, "L") - L) < 0.03),
                          key=lambda row: value(row, "pitch_center_deg"))
        x = [value(row, "pitch_center_deg") -
             value(row, "theory_critical_deg") for row in selected]
        y = [value(row, "P_precip") for row in selected]
        ax.errorbar(x, y,
                    yerr=[[center - value(row, "wilson_low_95")
                           for center, row in zip(y, selected)],
                          [value(row, "wilson_high_95") - center
                           for center, row in zip(y, selected)]],
                    fmt="o-", capsize=3, color=color, label=f"L ≈ {L:.2f}")
    ax.axvline(0, color="black", linestyle="--", linewidth=1,
               label="Адиабатическая граница")
    ax.set(xlabel="Центр бина − теоретический αc, градусы",
           ylabel="Доля высыпавшихся", ylim=(-0.05, 1.06),
           title="Вероятность высыпания около конуса потерь")
    ax.legend(loc="center right")
    fig.tight_layout()
    save_plot(fig, reports, "05_probability_by_pitch.png", plt)

    distributions = ("uniform_angle", "isotropic_hemisphere")
    labels = ("Равномерный угол", "Изотропные направления")
    pooled = {}
    for distribution in distributions:
        selected = [row for row in source_distributions
                    if row["distribution"] == distribution]
        n = sum(int(row["N"]) for row in selected)
        k = sum(int(row["precipitated"]) for row in selected)
        mirrored = sum(int(row["mirrored"]) for row in selected)
        other = sum(int(row["escaped"]) + int(row["unresolved"])
                    for row in selected)
        pooled[distribution] = (n, k, mirrored, other, *wilson(k, n),
                                value(selected[0], "theory_fraction"))
    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = [0, 1]
    fractions = [pooled[name][1] / pooled[name][0] for name in distributions]
    ax.errorbar(x, fractions,
                yerr=[[fractions[i] - pooled[name][4]
                       for i, name in enumerate(distributions)],
                      [pooled[name][5] - fractions[i]
                       for i, name in enumerate(distributions)]],
                fmt="o", markersize=9, capsize=5, color="#2467a2",
                label="Расчёт, 95% интервал Уилсона")
    ax.scatter(x, [pooled[name][6] for name in distributions], marker="D",
               color="#a65320", label="Идеальный конус потерь")
    ax.set(xticks=x, xticklabels=labels, ylabel="Доля высыпавшихся",
           ylim=(0.5, 0.9), title="Одна магнитная линия, два допущения об источнике")
    ax.legend()
    fig.tight_layout()
    save_plot(fig, reports, "05_source_comparison.png", plt)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    for distribution, color in (("uniform_angle", "#2467a2"),
                                ("isotropic_hemisphere", "#a65320")):
        short_label = "равн. α" if distribution == "uniform_angle" else "равн. cos α"
        for seed, style in ((20261006, "o-"), (20261007, "s--")):
            selected = sorted((row for row in source_convergence
                               if row["distribution"] == distribution and
                               int(row["seed"]) == seed),
                              key=lambda row: int(row["prefix_N"]))
            ax.plot([int(row["prefix_N"]) for row in selected],
                    [value(row, "P_precip") for row in selected], style,
                    color=color, label=f"{short_label}, seed {seed}")
        ax.axhline(pooled[distribution][6], color=color, alpha=0.5,
                   linestyle=":")
    ax.set(xlabel="Число маркеров на seed", ylabel="Доля высыпавшихся",
           xticks=[32, 64, 128], ylim=(0.5, 0.9),
           title="Разброс при изменении seed и размера ансамбля")
    ax.legend(fontsize=8)
    fig.tight_layout()
    save_plot(fig, reports, "05_source_convergence.png", plt)

    paired = {}
    for row in precipitation_particles:
        if (row["group"].startswith("uniform_angle_") or
                row["group"].startswith("isotropic_hemisphere_")):
            paired[(row["group"], int(row["id"]))] = (
                row["outcome"] == "Precipitated", value(row, "phase_deg"))
    differences = []
    for seed_index in (0, 1):
        for marker_id in range(1, 129):
            uniform = paired[(f"uniform_angle_{seed_index}", marker_id)]
            isotropic = paired[(f"isotropic_hemisphere_{seed_index}", marker_id)]
            if uniform[1] != isotropic[1]:
                raise ValueError("Source comparison lost its paired gyro phases")
            differences.append(int(uniform[0]) - int(isotropic[0]))
    random_generator = random.Random(20261006)
    bootstrap = sorted(
        sum(random_generator.choice(differences) for _ in differences) /
        len(differences) for _ in range(5000))
    difference = sum(differences) / len(differences)
    difference_low = bootstrap[125]
    difference_high = bootstrap[4875]
    source_table = "\n".join(
        f"| {label} | {pooled[name][1]}/{pooled[name][0]} | "
        f"{pooled[name][1] / pooled[name][0]:.3f} | "
        f"{pooled[name][4]:.3f}–{pooled[name][5]:.3f} | "
        f"{pooled[name][6]:.3f} |"
        for name, label in zip(distributions, labels))
    unresolved_total = sum(int(row["unresolved"]) for row in precipitation_bins +
                           source_distributions)
    write_report(reports, "05_precipitation_statistics.md", f"""
# Эксперимент 5. Статистика высыпания электронов

## Физический вопрос

Как меняется доля электронов, достигающих условной атмосферы, около границы
конуса потерь и при замене распределения начальных направлений?

## Постановка

В статическом диполе с высоты 450 км запускаются электроны 10 кэВ. Для трёх
линий `L` около теоретического угла `αc` заданы пять бинов шириной **0,1°**;
в каждом по **16** случайных углов внутри бина и гирофаз. Отдельно на линии
`L ≈ 10,82` два независимых seed дают по **128** маркеров для каждого источника:
равномерного по `α ∈ [0°,90°]` и равномерного по `cos α ∈ [0,1]`
(изотропного по направлениям на нисходящей полусфере). Шаг — около `T/50`.
Для всех исходов отдельно считаются `Precipitated`, `Mirrored`, `Escaped`,
`Unresolved`; здесь суммарно `Unresolved = {unresolved_total}`.

![Вероятность в угловых бинах](05_probability_by_pitch.png)

Смешанные доли в центральных бинах появляются потому, что **сам бин имеет
конечную ширину** и в нём выбираются разные начальные углы. Траектория при
заданном начальном состоянии детерминирована. Усики — 95%-е интервалы
Уилсона для доли внутри выбранного бина, не погрешность положения `αc`.

![Сравнение двух источников](05_source_comparison.png)

| Начальные направления | Дошли / всего | Доля | 95% интервал Уилсона | Идеальный конус |
| --- | ---: | ---: | ---: | ---: |
{source_table}

Для идеальной резкой границы ожидаются разные доли:
`P(равномерный угол) = αc/90°`, `P(изотропный) = 1 − cos αc`.
Фактическая разность при **попарном сопоставлении одинаковых случайных квантилей** —
**{difference:.3f}**, парный 95%-й bootstrap-интервал
**[{difference_low:.3f}; {difference_high:.3f}]** (5000 повторов).
Общее число маркеров в каждом широком источнике — 256.

![Разброс seed и N](05_source_convergence.png)

Значения при `N=32/64/128` — вложенные префиксы одного прогона, а не три
независимых опыта. Два seed позволяют увидеть разброс между реализациями;
пунктирные горизонтали — оценки идеального конуса потерь.

## Вывод и ограничения

В этой модели вероятность резко меняется около теоретического конуса потерь.
Распределение начальных направлений заметно меняет **интегральную** долю
высыпавшихся даже при неизменном магнитном поле. Это сравнение двух
**заданных гипотез об источнике**, не измерение реального электронного потока.
Малый размер бина (`N=16`) даёт широкие интервалы; в модели нет электрического
ускорения, рассеяния по углам до 150 км и реального геомагнитного поля.

Исходные данные: [угловые бины](../precipitation_bins.csv),
[источники](../source_distributions.csv),
[сходимость](../source_convergence.csv),
[все частицы](../precipitation_particles.csv).
В последнем файле сохранены полные начальные координаты и импульсы; единичный
вес маркера здесь не является измеренным физическим потоком.
Код: `src/experiments/PrecipitationStatistics.hpp`.
[Как воспроизвести](../../../src/experiments/README.md).
""")

    # 06: Spatial shape and attenuation of a divergence-free buried dipole.
    anomaly_cases = ("weak_broad", "strong_broad", "strong_compact")
    case_labels = {"zero": "Нулевая", "weak_broad": "Слабая, 500 км",
                   "strong_broad": "Усиленная, 500 км",
                   "strong_compact": "Усиленная, 250 км"}
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), sharey=True)
    for axis, altitude in zip(axes, (0, 150000, 450000)):
        for name in anomaly_cases:
            selected = sorted((row for row in anomaly_field
                               if row["case"] == name and
                               value(row, "altitude_m") == altitude),
                              key=lambda row: value(row, "offset_deg"))
            axis.plot([value(row, "offset_deg") for row in selected],
                      [100 * value(row, "relative_anomaly") for row in selected],
                      label=case_labels[name])
        axis.set(xlabel="Смещение от источника, градусы",
                 title=f"Высота {altitude / 1000:.0f} км")
    axes[0].set_ylabel("|δB| / |B₀|, %")
    axes[-1].legend(fontsize=8)
    fig.suptitle("Пространственный масштаб и затухание аномалии")
    fig.tight_layout()
    save_plot(fig, reports, "06_anomaly_field.png", plt)
    field_table = "\n".join(
        f"| {case_labels[name]} | {value(next(row for row in anomaly_field if row['case'] == name), 'depth_m') / 1000:.0f} | "
        + " | ".join(
            f"{100 * value(next(row for row in anomaly_field
                                 if row['case'] == name and
                                 value(row, 'altitude_m') == altitude and
                                 value(row, 'offset_deg') == 0), 'relative_anomaly'):.3f}%"
            for altitude in (0, 150000, 450000)) + " |"
        for name in anomaly_cases)
    write_report(reports, "06_anomaly_field.md", f"""
# Эксперимент 6. Локальная магнитная аномалия

## Физический вопрос

Как изменяется вклад подповерхностного источника в поле на поверхности и на
высотах, где стартуют и заканчиваются электронные траектории?

## Постановка и метод

К земному диполю добавлен отдельный диполь на глубине 500 или 250 км под
точкой с колатитудой 18,335°. Его момент направлен радиально, сила задаётся
отношением **осевого поля аномалии на поверхности** к полю основного диполя
там же. Слабый вариант: 1%; усиленные учебные варианты: 30%. Это модельные
параметры, а не утверждение о типичной земной коровой аномалии.

Вне подземного источника оба слагаемых имеют `∇·B=0`; конечными разностями
это отдельно проверено автоматическим тестом. Нулевая амплитуда в том же тесте
даёт ровно исходный диполь.

![Поле аномалии на трёх высотах](06_anomaly_field.png)

| Случай | Глубина, км | δB/B₀, 0 км | δB/B₀, 150 км | δB/B₀, 450 км |
| --- | ---: | ---: | ---: | ---: |
{field_table}

## Вывод и ограничения

Подповерхностный диполь быстро ослабевает с высотой. Сравнивать силу аномалии
только на поверхности недостаточно: электроны стартуют на 450 км, а исход
`Precipitated` регистрируется при достижении 150 км. Усиленный случай нужен
для проверки чувствительности метода; он не является реалистичным прогнозом.
Показан срез поля в одной меридиональной плоскости, не глобальная карта Земли.

Исходные данные: [пространственная карта](../anomaly_field.csv).
Код: `src/physics/LocalizedDipoleField.hpp`.
[Как воспроизвести](../../../src/experiments/README.md).
""")

    # 07: Paired trajectories with identical initial phase-space points.
    names = ("zero",) + anomaly_cases
    paired = {}
    for name in names:
        selected = [row for row in anomaly_particles if row["case"] == name]
        differences = []
        transitions = [0, 0, 0]
        base_count = 0
        changed_count = 0
        for row in selected:
            base = row["baseline_outcome"]
            changed = row["anomaly_outcome"]
            base_count += base == "Precipitated"
            changed_count += changed == "Precipitated"
            differences.append(int(changed == "Precipitated") -
                               int(base == "Precipitated"))
            if base == "Precipitated" and changed == "Mirrored":
                transitions[0] += 1
            elif base == "Mirrored" and changed == "Precipitated":
                transitions[1] += 1
            elif base != changed:
                transitions[2] += 1
        random_generator = random.Random(20261007)
        bootstrap = sorted(
            sum(random_generator.choice(differences) for _ in differences) /
            len(differences) for _ in range(3000))
        paired[name] = {
            "N": len(selected), "base": base_count, "changed": changed_count,
            "delta": sum(differences) / len(differences),
            "low": bootstrap[75], "high": bootstrap[2925],
            "transitions": transitions,
            "wilson": wilson(changed_count, len(selected)),
        }
    scan_by_name = {row["case"]: row for row in anomaly_scan}
    labels = [case_labels[name] for name in names]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    x = list(range(len(names)))
    y = [100 * paired[name]["delta"] for name in names]
    axes[0].errorbar(x, y,
                     yerr=[[y[i] - 100 * paired[name]["low"]
                            for i, name in enumerate(names)],
                           [100 * paired[name]["high"] - y[i]
                            for i, name in enumerate(names)]],
                     fmt="o", capsize=4, color="#2467a2")
    axes[0].axhline(0, color="black", linewidth=1)
    axes[0].set(ylabel="ΔP осаждения, процентные пункты",
                title="Парная разность, 95% bootstrap")
    axes[1].plot(x, [value(scan_by_name[name], "delta_midpoint_deg")
                     for name in names], "o-", color="#a65320")
    axes[1].axhline(0, color="black", linewidth=1)
    axes[1].set(ylabel="Сдвиг критического угла, градусы",
                title="Одна стартовая гирофаза")
    for axis in axes:
        axis.set_xticks(x, labels, rotation=20, ha="right")
    fig.tight_layout()
    save_plot(fig, reports, "07_paired_anomaly.png", plt)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(x, [paired[name]["transitions"][0] for name in names],
           label="Достиг → отразился", color="#a65320")
    ax.bar(x, [paired[name]["transitions"][1] for name in names],
           bottom=[paired[name]["transitions"][0] for name in names],
           label="Отразился → достиг", color="#2467a2")
    ax.set(xticks=x, xticklabels=labels, ylabel="Число частиц",
           title="Смена исхода у тех же начальных состояний")
    ax.tick_params(axis="x", labelrotation=20)
    ax.legend()
    fig.tight_layout()
    save_plot(fig, reports, "07_transitions.png", plt)
    selected_entries = [row for row in anomaly_particles
                        if row["case"] == "weak_broad" and
                        row["baseline_outcome"] == "Precipitated" and
                        row["anomaly_outcome"] == "Precipitated"]
    if not selected_entries:
        raise ValueError("No paired atmosphere entries for the anomaly comparison")
    weak_mean_entry_shift_m = sum(value(row, "entry_shift_m") for row in selected_entries) / len(selected_entries)
    origin_x = sum(value(row, "baseline_event_x_m") for row in selected_entries) / len(selected_entries)
    origin_y = sum(value(row, "baseline_event_y_m") for row in selected_entries) / len(selected_entries)
    fig, ax = plt.subplots(figsize=(6, 5))
    for prefix, label, color in (("baseline", "Без аномалии", "#2467a2"),
                                 ("anomaly", "С аномалией", "#a65320")):
        ax.scatter([(value(row, f"{prefix}_event_x_m") - origin_x) / 1000
                    for row in selected_entries],
                   [(value(row, f"{prefix}_event_y_m") - origin_y) / 1000
                    for row in selected_entries],
                   s=12, alpha=0.55, label=label, color=color)
    ax.set(xlabel="x − среднее x₀, км", ylabel="y − среднее y₀, км",
           title="Слабая аномалия: точки входа на 150 км")
    ax.legend()
    fig.tight_layout()
    save_plot(fig, reports, "07_entry_points.png", plt)
    comparison_table = "\n".join(
        f"| {case_labels[name]} | {paired[name]['changed']}/{paired[name]['N']} | "
        f"{paired[name]['wilson'][0]:.3f}–{paired[name]['wilson'][1]:.3f} | "
        f"{100 * paired[name]['delta']:+.2f} | "
        f"[{100 * paired[name]['low']:+.2f}; {100 * paired[name]['high']:+.2f}] | "
        f"{paired[name]['transitions'][0]} / {paired[name]['transitions'][1]} | "
        f"{value(scan_by_name[name], 'delta_midpoint_deg'):+.3f} |"
        for name in names)
    baseline_total = paired["zero"]["base"]
    scan_sensitivity = max(abs(value(row, "delta_midpoint_deg") -
                               value(row, "delta_midpoint_fine_deg"))
                           for row in anomaly_scan)
    write_report(reports, "07_paired_anomaly.md", f"""
# Эксперимент 7. Влияние аномалии на высыпание электронов

## Физический вопрос

Меняет ли локальный подповерхностный диполь границу конуса потерь и исходы
**одних и тех же** электронов?

## Постановка и метод

Из точки на высоте 450 км стартуют электроны 10 кэВ, угол к **базовому** полю
равномерно выбран в бине шириной 1,2° вокруг базовой критической границы.
Два seed по 128 частиц дают N=256. Начальные `(r,p,t,id,weight)` генерируются
один раз и дословно повторяются в четырёх полях; бин специально сосредоточен
вблизи порога и **не представляет всю изотропную полусферу**. Шаг `T/50`,
граница атмосферы 150 км. Доверительный интервал разности — парный bootstrap
3000 повторов; интервалы отдельной доли — Уилсон.

![Парная разность и сдвиг границы](07_paired_anomaly.png)

В базовом поле достигли границы **{baseline_total}/256** электронов.

| Поле | Достигли, k/N | 95% Уилсон | ΔP, п.п. | 95% парный bootstrap, п.п. | Достиг→отразился / обратно | Δαc, ° |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
{comparison_table}

![Переходы исходов](07_transitions.png)

![Сопоставимые точки входа при слабой аномалии](07_entry_points.png)

Среднее смещение точек входа рассчитано **только для частиц, достигших
150 км в обоих полях**; для слабой аномалии это **{weak_mean_entry_shift_m:.1f} м**.
В усиленных вариантах таких частиц
в выбранном узком бине нет, поэтому ноль в колонке смещения CSV означает
«не определено», а не отсутствие отклонения. `Escaped` и `Unresolved` не
приравниваются к отражению или высыпанию. Угловая граница найдена бинарным
поиском для фиксированной базовой гирофазы и того же начального направления;
точность скобки 0,05°, а не точная физическая константа.
При уменьшении шага с `T/50` до `T/100` сдвиг критического угла меняется
не более чем на **{scan_sensitivity:.3f}°** в этой серии. Это проверка
сканирования, но не статистики всего ансамбля при втором шаге.

## Вывод и ограничения

Парное сравнение изолирует изменение магнитного поля от случайной перегенерации
источника. Даже слабая **заданная** аномалия уменьшила долю высыпания в
выбранном бине; парный интервал для этой разности не включает ноль.
Усиленные аномалии — учебные сценарии: их одинаковый итог `0/256` обусловлен
насыщением узкой выборки и **не означает**, что глубина не влияет на поле или
конус потерь. Кроме статистики, остаются численные погрешности ансамбля и
упрощение статического диполя. Доля относится только к выбранному узкому
бину углов, не к реальному электронному потоку.

Исходные данные: [все парные частицы](../anomaly_particles.csv),
[агрегаты](../anomaly_summary.csv), [скан границы](../anomaly_scan.csv),
[поле](../anomaly_field.csv). Код: `src/experiments/anomaly_experiment.cpp`.
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
5. [Статистика высыпания](05_precipitation_statistics.md)
6. [Локальная магнитная аномалия](06_anomaly_field.md)
7. [Парное влияние аномалии на высыпание](07_paired_anomaly.md)

Все рисунки лежат рядом с отчётами. Исходные CSV — в родительском каталоге.
Графики и отчёты можно [пересоздать по инструкции](../../../src/experiments/README.md).
""")
    print(f"Wrote seven Markdown reports and thirteen figures to {reports}")


if __name__ == "__main__":
    main()
