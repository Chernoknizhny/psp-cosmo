#!/usr/bin/env python3
"""
Построение Рис. 1 для Документа №05 (ред. 4):
Распределение доли газа в гало для изолированных и взаимодействующих галактик.

Данные — из оригинала PDF, Таблица 13 (45 измерений).
Взаимодействующих: 13. Изолированных: 32.

Использование: python3 make_plot.py
Выход: combined_halo_plot_histogram.png, data_halo_gas.csv
"""

import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

data = [
    (1,  "NGC 4449",   59.4, "Взаимодействующая"),
    (2,  "NGC 4532",   54.4, "Взаимодействующая"),
    (3,  "NGC 3169",   51.8, "Взаимодействующая"),
    (4,  "NGC 3368",   42.2, "Взаимодействующая"),
    (5,  "NGC 5194",   36.5, "Взаимодействующая"),
    (6,  "NGC 5194*",  36.5, "Взаимодействующая"),
    (7,  "NGC 4490",   34.7, "Взаимодействующая"),
    (8,  "IC 1727",    25.3, "Взаимодействующая"),
    (9,  "NGC 4254",   24.5, "Взаимодействующая"),
    (10, "NGC 4725",   23.8, "Взаимодействующая"),
    (11, "NGC 5033",   20.6, "Взаимодействующая"),
    (12, "NGC 672",    17.7, "Взаимодействующая"),
    (13, "NGC 5457*",  12.7, "Взаимодействующая"),
    (14, "NGC 2841",   40.1, "Изолированная"),
    (15, "NGC 2841*",  40.0, "Изолированная"),
    (16, "NGC 4414",   30.3, "Изолированная"),
    (17, "NGC 5055",   28.1, "Изолированная"),
    (18, "NGC 5055*",  28.0, "Изолированная"),
    (19, "NGC 3718",   25.0, "Изолированная"),
    (20, "NGC 7331*",  25.0, "Изолированная"),
    (21, "NGC 4618",   23.7, "Изолированная"),
    (22, "NGC 3338",   22.0, "Изолированная"),
    (23, "NGC 4496",   21.5, "Изолированная"),
    (24, "NGC 5248",   20.3, "Изолированная"),
    (25, "NGC 3344",   19.3, "Изолированная"),
    (26, "NGC 2903",   19.1, "Изолированная"),
    (27, "NGC 2903*",  19.0, "Изолированная"),
    (28, "NGC 925*",   18.0, "Изолированная"),
    (29, "NGC 2541",   17.0, "Изолированная"),
    (30, "NGC 4536",   16.9, "Изолированная"),
    (31, "NGC 3486",   16.7, "Изолированная"),
    (32, "NGC 4303",   16.7, "Изолированная"),
    (33, "NGC 5474",   16.5, "Изолированная"),
    (34, "Sex B",      16.4, "Изолированная"),
    (35, "NGC 4214",   16.0, "Изолированная"),
    (36, "NGC 3198*",  16.0, "Изолированная"),
    (37, "NGC 864",    15.9, "Изолированная"),
    (38, "NGC 3198",   15.8, "Изолированная"),
    (39, "NGC 628*",   14.0, "Изолированная"),
    (40, "NGC 628",    13.7, "Изолированная"),
    (41, "NGC 3521*",  13.0, "Изолированная"),
    (42, "NGC 3521",   12.6, "Изолированная"),
    (43, "NGC 4395",    9.9, "Изолированная"),
    (44, "NGC 5457",    9.3, "Изолированная"),
    (45, "NGC 4559",    6.8, "Изолированная"),
]

iso = [d[2] for d in data if d[3] == "Изолированная"]
intg = [d[2] for d in data if d[3] == "Взаимодействующая"]

# --- Статистика ---
def stats(arr, label):
    a = np.array(arr)
    print(f"{label} (N={len(a)}): среднее {np.mean(a):.1f}%, "
          f"медиана {np.median(a):.1f}%, ст.откл {np.std(a, ddof=1):.1f}%, "
          f"мин {np.min(a):.1f}%, макс {np.max(a):.1f}%")

stats(iso,  "Изолированные")
stats(intg, "Взаимодействующие")

# --- CSV ---
with open("data_halo_gas.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["N", "Galaxy", "f_halo_%", "Type"])
    for row in data:
        w.writerow(row)

# --- График ---
plt.rcParams.update({"font.size": 11, "figure.dpi": 150})

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.5),
                                gridspec_kw={"width_ratios": [1.4, 1]})

c_iso, c_intg = "#3274C1", "#E1812C"

# Гистограмма
bins = np.arange(0, 70, 5)
ax1.hist([iso, intg], bins=bins, color=[c_iso, c_intg],
         label=[f"Изолированные (N={len(iso)})",
                f"Взаимодействующие (N={len(intg)})"],
         edgecolor="white", linewidth=0.6, alpha=0.85)
ax1.set_xlabel("Доля газа в гало, %")
ax1.set_ylabel("Число галактик")
ax1.set_title("Распределение доли газа в гало")
ax1.legend(loc="upper right", fontsize=10)
ax1.set_xlim(0, 65)
ax1.grid(axis="y", alpha=0.3, linestyle="--")
ax1.set_axisbelow(True)

# Box-plot
bp = ax2.boxplot([iso, intg], tick_labels=["Изолир.", "Взаим."],
                 patch_artist=True, widths=0.5,
                 medianprops=dict(color="black", linewidth=2))
for patch, color in zip(bp["boxes"], [c_iso, c_intg]):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
    patch.set_edgecolor("black")

for i, (arr, color) in enumerate(zip([iso, intg], [c_iso, c_intg]), 1):
    jitter = np.random.default_rng(42).uniform(-0.12, 0.12, len(arr))
    ax2.scatter(np.full(len(arr), i) + jitter, arr,
                color=color, alpha=0.5, s=20, zorder=3,
                edgecolors="white", linewidths=0.3)
ax2.set_ylabel("Доля газа в гало, %")
ax2.set_title("Диаграмма размаха")
ax2.grid(axis="y", alpha=0.3, linestyle="--")
ax2.set_axisbelow(True)

fig.suptitle("Рис. 1. Газовое гало галактик: изолированные vs взаимодействующие\n"
             "(объединённая выборка, N=45; данные Wang et al. 2025 + FEASTS)",
             fontsize=13, y=0.98)
plt.tight_layout(rect=[0, 0, 1, 0.93])
plt.savefig("combined_halo_plot_histogram.png", dpi=150,
            bbox_inches="tight", facecolor="white")
print("Готово: combined_halo_plot_histogram.png")
