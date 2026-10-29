"""Simulate Turbie for the provided wind cases and save the results and plots."""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

import turbie_mod as tm


# --- Input and output folders ---
wind_folder = Path("inputs/wind_files")
output_folder = Path("outputs")
output_folder.mkdir(exist_ok=True)

turbie_parameters_file = "inputs/turbie_inputs/turbie_parameters.txt"
C_T_table_file = "inputs/turbie_inputs/CT.txt"

# --- Fixed setup (done once) ---
# Turbie parameters and matrices
parameters = tm.read_turbie_parameters(turbie_parameters_file)
mb, mn, mh, mt, c1, c2, k1, k2, fb, ft, drb, drt, Dr, rho = tm.get_turbie_parameters(parameters)
M, C, K = tm.calculate_MCK_matrix(mb, mn, mh, mt, c1, c2, k1, k2)

# Rotor area and thrust coefficient table
A = tm.calculate_rotor_area(Dr)
C_T_table = tm.read_C_T_table(C_T_table_file)

# Initial state [x1, x2, x1_dot, x2_dot]: Turbie starts at rest
y0 = [0, 0, 0, 0]

# Initial transient discarded when computing statistics [s]
t_transient = 60.0

# Statistics of all TI categories, used for the plots
all_statistics = {}

# --- Loop over TI categories ---
for ti_folder in sorted(wind_folder.iterdir()):
    print("TI category:", ti_folder.name)

    # Wind files of this TI category, sorted by wind speed
    wind_files = sorted(ti_folder.glob("*.txt"), key=lambda f: float(f.stem.split("_")[1]))

    # Statistics of this TI category: one row per wind speed
    statistics = []

    # --- Loop over wind speed files ---
    for wind_file in wind_files:
        print("  Simulating", wind_file.name)

        # Wind time series and thrust coefficient
        time, wind_speed = tm.read_wind_file(wind_file)
        mean_wind_speed = np.mean(wind_speed)
        C_T = tm.get_C_T(mean_wind_speed, C_T_table)

        # Solve the dynamics of the system
        sol = solve_ivp(
            tm.state_derivative,
            t_span=[time[0], time[-1]],
            y0=y0,
            t_eval=time,
            args=(M, C, K, time, wind_speed, C_T, rho, A),
        )

        # Save the time series of this case
        results_file = output_folder / ti_folder.name / f"results_{wind_file.stem}.txt"
        results = np.column_stack((sol.t, wind_speed, sol.y[0], sol.y[1], sol.y[2], sol.y[3]))
        np.savetxt(
            results_file,
            results[::10],  # Save every 10th point to reduce file size
            fmt=["%.2f", "%.3f", "%.5f", "%.5f", "%.5f", "%.5f"],
            header="time [s], wind speed [m/s], x1 [m], x2 [m], x1_dot [m/s], x2_dot [m/s]",
        )

        # Statistics without the initial transient
        mask = sol.t >= t_transient
        mean_x1, std_x1 = np.mean(sol.y[0][mask]), np.std(sol.y[0][mask])
        mean_x2, std_x2 = np.mean(sol.y[1][mask]), np.std(sol.y[1][mask])
        statistics.append([mean_wind_speed, mean_x1, std_x1, mean_x2, std_x2])

        # Time series plot for one representative case
        if wind_file.stem == "wind_6_ms_TI_0.1":
            fig, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
            fig.suptitle("Turbie response, mean wind speed 6 m/s, TI = 0.1")

            # Panel 1: wind speed vs time
            axs[0].plot(time, wind_speed)
            axs[0].set_ylabel("Wind Speed [m/s]")

            # Panel 2: blade displacement x1 vs time
            axs[1].plot(time, sol.y[0])
            axs[1].set_ylabel("Blade Displacement [m]")

            # Panel 3: tower displacement x2 vs time
            axs[2].plot(time, sol.y[1])
            axs[2].set_ylabel("Tower Displacement [m]")
            axs[2].set_xlabel("Time [s]")
            axs[2].set_xlim(time[0], time[-1])

            # Common formatting
            for ax in axs:
                ax.grid(True)
            fig.tight_layout()

            # Save the figure and close it
            fig.savefig(output_folder / "time_series_6_ms_TI_0.1.png", dpi=300)
            plt.close(fig)

    # Save the statistics of this TI category
    statistics = np.array(statistics)
    statistics_file = output_folder / ti_folder.name / f"statistics_{ti_folder.name}.txt"
    np.savetxt(
        statistics_file,
        statistics,
        fmt="%.6f",
        header="mean wind speed [m/s], mean x1 [m], std x1 [m], mean x2 [m], std x2 [m]",
    )

    # Keep the statistics of this TI category for the plots
    all_statistics[ti_folder.name] = statistics


# --- Statistics plot: mean and standard deviation vs mean wind speed ---
fig, axs = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
fig.suptitle("Turbie displacement statistics vs mean wind speed")

# One curve per TI category in each panel
for ti_name, statistics in all_statistics.items():
    ti_label = "TI = " + ti_name.split("_")[-1]
    wind = statistics[:, 0]

    axs[0, 0].plot(wind, statistics[:, 1], marker="o", label=ti_label)  # mean x1
    axs[1, 0].plot(wind, statistics[:, 2], marker="o", label=ti_label)  # std x1
    axs[0, 1].plot(wind, statistics[:, 3], marker="o", label=ti_label)  # mean x2
    axs[1, 1].plot(wind, statistics[:, 4], marker="o", label=ti_label)  # std x2

# Titles of the columns and labels of the axes
axs[0, 0].set_title("Blades (x1)")
axs[0, 1].set_title("Tower (x2)")
axs[0, 0].set_ylabel("Mean displacement [m]")
axs[1, 0].set_ylabel("Std of displacement [m]")
axs[1, 0].set_xlabel("Mean wind speed [m/s]")
axs[1, 1].set_xlabel("Mean wind speed [m/s]")

for ax in axs.flat:
    ax.grid(True)
axs[0, 0].legend()
fig.tight_layout()

# Save the figure and close it
fig.savefig(output_folder / "statistics_vs_wind_speed.png", dpi=300)
plt.close(fig)