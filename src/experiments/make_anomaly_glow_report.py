#!/usr/bin/env python3
"""Illustrated paired anomaly -> precipitation -> emission report."""

import atexit
import csv
import math
import os
from pathlib import Path
import random
import sys
import tempfile


def load(path):
    with path.open(newline="", encoding="utf-8") as stream:
        data = list(csv.DictReader(stream))
    if not data:
        raise ValueError(f"No data in {path}")
    return data


def paired_interval(deltas, seed=6413, repetitions=10000):
    rng = random.Random(seed)
    count = len(deltas)
    distribution = sorted(sum(deltas[rng.randrange(count)] for _ in range(count)) /
                          count for _ in range(repetitions))
    return distribution[int(0.025 * repetitions)], distribution[int(0.975 * repetitions)]


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
        raise SystemExit("Usage: python3 src/experiments/make_anomaly_glow_report.py RESULT_DIRECTORY")
    directory = Path(sys.argv[1]).resolve()
    reports = directory / "reports"
    reports.mkdir(parents=True, exist_ok=True)
    particles = load(directory / "anomaly_glow_particles.csv")
    summary = load(directory / "anomaly_glow_summary.csv")
    profile = load(directory / "anomaly_glow_profile.csv")
    totals = load(directory / "anomaly_glow_totals.csv")
    diagnostics = load(directory / "anomaly_glow_diagnostics.csv")
    seeds = sorted({int(row["seed"]) for row in particles})
    if len(seeds) < 2:
        raise ValueError("At least two independent seeds are required")
    if any(row["zero_50"] != row["baseline_50"] for row in particles):
        raise ValueError("Zero-anomaly outcome control failed")
    if any(row["baseline_50"] != row["baseline_100"] or
           row["weak_50"] != row["weak_100"] for row in particles):
        raise ValueError("Outcome changed after dt was halved")
    if any(int(row["escaped"]) or int(row["unresolved"]) for row in summary):
        raise ValueError("Escaped or unresolved particle: do not treat as mirrored")

    baseline = [row["baseline_50"] == "Precipitated" for row in particles]
    weak = [row["weak_50"] == "Precipitated" for row in particles]
    deltas = [int(after) - int(before) for before, after in zip(baseline, weak)]
    n = len(particles)
    base_count = sum(baseline)
    weak_count = sum(weak)
    lost = sum(before and not after for before, after in zip(baseline, weak))
    gained = sum(after and not before for before, after in zip(baseline, weak))
    low, high = paired_interval(deltas)
    base_low, base_high = wilson(base_count, n)
    weak_low, weak_high = wilson(weak_count, n)
    by_seed = {seed: [row for row in particles if int(row["seed"]) == seed]
               for seed in seeds}
    for row in summary:
        seed = int(row["seed"])
        if seed not in by_seed:
            raise ValueError(f"Unknown seed in summary: {seed}")
        scenario = row["scenario"]
        step = row["steps_per_gyroperiod"]
        column = {("baseline", "50"): "baseline_50",
                  ("zero", "50"): "zero_50",
                  ("weak", "50"): "weak_50",
                  ("baseline", "100"): "baseline_100",
                  ("weak", "100"): "weak_100"}.get((scenario, step))
        if column is None:
            raise ValueError(f"Unknown scenario/step: {scenario}/{step}")
        expected = sum(item[column] == "Precipitated" for item in by_seed[seed])
        if expected != int(row["precipitated"]) or len(by_seed[seed]) != int(row["N"]):
            raise ValueError(f"Summary disagrees with paired particles: {seed}/{column}")

    cache = tempfile.TemporaryDirectory(prefix="aurora-matplotlib-")
    atexit.register(cache.cleanup)
    os.environ["MPLCONFIGDIR"] = cache.name
    os.environ["XDG_CACHE_HOME"] = cache.name
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.grid": True, "grid.alpha": 0.25})
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    fractions = []
    for seed in seeds:
        batch = by_seed[seed]
        fractions.append((sum(row["baseline_50"] == "Precipitated" for row in batch) /
                          len(batch),
                          sum(row["weak_50"] == "Precipitated" for row in batch) /
                          len(batch)))
    for index, (first, second) in enumerate(fractions):
        axes[0].plot([0, 1], [first, second], "o-", alpha=0.7,
                     label=f"seed {seeds[index]}")
    axes[0].set_xticks([0, 1], ["Диполь", "Слабая аномалия"])
    axes[0].set(ylabel="Доля достигших 150 км", ylim=(0, 1),
                title="Парные доли по seed")
    axes[0].legend(fontsize=8)
    per_seed_delta = [(second - first) * 100 for first, second in fractions]
    axes[1].bar(range(len(seeds)), per_seed_delta, color="#7b3e92")
    axes[1].axhline(0, color="black", linewidth=0.8)
    axes[1].set_xticks(range(len(seeds)), [str(seed) for seed in seeds], rotation=20)
    axes[1].set(ylabel="ΔP, процентные пункты", title="Слабая аномалия − диполь")
    fig.tight_layout()
    fig.savefig(reports / "08_anomaly_precipitation.png", dpi=170)
    plt.close(fig)

    color_by_line = {"O_5577": "#258243", "O_6300": "#b94d43"}
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
    for axis, line in zip(axes, color_by_line):
        for scenario, style in (("baseline", "-"), ("weak", "--")):
            per_height = {}
            for row in profile:
                if row["scenario"] != scenario or row["line"] != line:
                    continue
                height = int(row["height_bin_km"])
                per_height[height] = per_height.get(height, 0.0) + float(row["photon_rate_s1"])
            if not per_height:
                raise ValueError(f"No profile for {scenario}, {line}")
            heights = sorted(per_height)
            label = "Базовое поле" if scenario == "baseline" else "Слабая аномалия"
            axis.plot(heights, [per_height[height] / len(seeds) for height in heights],
                      style, color=color_by_line[line], linewidth=2, label=label)
        wavelength = {"O_5577": "557,7 нм", "O_6300": "630,0 нм"}[line]
        axis.set(xlabel="Высота, км", title=f"{line} ({wavelength})")
        axis.legend()
    axes[0].set_ylabel("Фотонов/с в высотном слое на 1 e⁻/с")
    fig.suptitle("Свечение: среднее по трём seed; одинаковая сетка в каждой паре")
    fig.tight_layout()
    fig.savefig(reports / "08_anomaly_emission_profile.png", dpi=170)
    plt.close(fig)

    aggregate = {}
    profile_totals = {}
    for row in profile:
        key = int(row["seed"]), row["scenario"], row["line"]
        profile_totals[key] = profile_totals.get(key, 0.0) + float(row["photon_rate_s1"])
    for row in totals:
        key = row["scenario"], row["line"]
        aggregate.setdefault(key, []).append(float(row["photon_rate_s1"]))
        full_key = int(row["seed"]), row["scenario"], row["line"]
        if not math.isclose(profile_totals.get(full_key, 0.0),
                            float(row["photon_rate_s1"]), rel_tol=1e-10):
            raise ValueError(f"Profile does not sum to total for {full_key}")
    for key, values in aggregate.items():
        if len(values) != len(seeds):
            raise ValueError(f"Incomplete emission totals for {key}")
    table = "\n".join(
        f"| {line} | {sum(aggregate[('baseline', line)]) / len(seeds):.3g} | "
        f"{sum(aggregate[('weak', line)]) / len(seeds):.3g} | "
        f"{(sum(aggregate[('weak', line)]) / sum(aggregate[('baseline', line)]) - 1) * 100:+.1f}% |"
        for line in color_by_line)
    changed_step = sum(row["baseline_50"] != row["baseline_100"] or
                       row["weak_50"] != row["weak_100"] for row in particles)
    if len(diagnostics) != 2 * len(seeds):
        raise ValueError("Incomplete atmospheric diagnostics")
    max_energy_residual = max(abs(float(row["energy_residual_eV_s1"])) /
                              float(row["initial_energy_eV_s1"])
                              for row in diagnostics)
    max_outside = max(float(row["outside_grid_rate_s1"]) for row in diagnostics)
    report = f"""# Эксперимент 8. От магнитной аномалии к синтетическому свечению

## Физический вопрос

Меняет ли слабый подповерхностный магнитный диполь высыпание **тех же** электронов,
и как это отражается на расчётном свечении кислорода?

## Постановка

Из одной точки на высоте 450 км стартуют электроны 10 кэВ: {len(seeds)} независимых
seed × {n // len(seeds)} частиц. Угол к **базовому** полю равномерно распределён в
узком бине 69,0 ± 0,6°. Начальные состояния и веса повторяются без изменений.
Контроль — нулевая аномалия; тест — локальный диполь на глубине 500 км с
максимальным полем на поверхности 1% от фонового поля в точке над источником.
Это *заданный модельный параметр*, а не измеренное поле конкретного места.
Сначала траектории доводятся до границы атмосферы 150 км или первого зеркала;
достигшие границы передаются атмосферному каскаду. Для парного расчёта свечения
используются общая объёмная сетка, одинаковые таблицы, один seed каскада и
суммарный вес источника 1 электрон/с. Графики показывают среднее по seed.

![Парные исходы](08_anomaly_precipitation.png)

Из {n} электронов в базовом поле достигли границы **{base_count}**,
при слабой аномалии — **{weak_count}**. Переходы: достиг→отразился **{lost}**,
отразился→достиг **{gained}**. Парная разность долей —
**{sum(deltas) / n * 100:+.1f} п.п.**; 95% парный bootstrap по частицам:
**[{low * 100:+.1f}; {high * 100:+.1f}] п.п.**. Индивидуальные разности по seed:
{', '.join(f'{item:+.1f}' for item in per_seed_delta)} п.п.
Отдельные 95% интервалы Уилсона: базовое поле
**[{base_low * 100:.1f}; {base_high * 100:.1f}]%**, слабая аномалия
**[{weak_low * 100:.1f}; {weak_high * 100:.1f}]%**.

Нулевой контроль совпал с базовым полем для всех исходов и конечных состояний.
При уменьшении шага Бориса с T/50 до T/100 классификация изменилась у
**{changed_step}/{n}** частиц. Это проверяет устойчивость *исходов* в выбранной
выборке; точность каждой координаты, каскада и интегральной яркости отдельно
не доказана.

## Высотный профиль свечения

![Высотные профили линий кислорода](08_anomaly_emission_profile.png)

Интеграл фотонной скорости по объёму, в фотонах/с при потоке источника 1 e⁻/с:

| Линия | Базовое поле | Слабая аномалия | Изменение |
| --- | ---: | ---: | ---: |
{table}

Это расчёт **эмиссии в объёме**, не яркости камеры и не видимой зелени GIF.
Изменение зависит и от числа проникших электронов, и от их точек входа и
направлений; три seed демонстрируют разброс, но не дают отдельного надёжного
интервала для фотонной разности. Парный интервал выше относится **только к
вероятности высыпания**.
Максимальный относительный остаток энергетического баланса каскадов в шести
запусках — **{max_energy_residual:.3g}**; потеря учтённых возбуждений за
границей общей сетки — **{max_outside:.3g} событий/с**. Высотные профили
дискретны по ячейкам сетки; линия между точками на графике — только для чтения.

## Вывод и пределы достоверности

В этой узкой группе около границы конуса потерь слабая модельная аномалия
меняет исходы и рассчитанное свечение; нулевой контроль и удвоение частоты
шага не указывают на артефакт классификации. Результат **не относится ко всему
электронному потоку**: углы специально сконцентрированы около порога, атмосфера
взята из одного профиля, а таблицы сечений имеют модельное продолжение до
25 кэВ. Модель даёт синтетическое, а не количественно предсказанное реальное
полярное сияние. GIF отдельно визуализирует траектории и постепенно раскрывает
уже рассчитанную эмиссию; он **не моделирует временную кинетику** заселённостей
и фотонов.

Исходные данные: [парные частицы](../anomaly_glow_particles.csv),
[сводка по шагу и seed](../anomaly_glow_summary.csv),
[высотный профиль](../anomaly_glow_profile.csv),
[интегральные скорости](../anomaly_glow_totals.csv).
[Диагностика каскадов](../anomaly_glow_diagnostics.csv),
[параметры прогона](../anomaly_glow_manifest.txt).
Разреженные 3D-поля сохранены рядом как `anomaly_glow_<seed>_<scenario>.csv`.
Код: `src/experiments/anomaly_glow_experiment.cpp`.
[Как воспроизвести](../../../src/experiments/README.md).
"""
    (reports / "08_anomaly_glow.md").write_text(report, encoding="utf-8")
    print(f"Saved {reports / '08_anomaly_glow.md'}")


if __name__ == "__main__":
    main()
