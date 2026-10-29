"""Simulate Turbie for the provided wind cases and save the results and plots."""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

import turbie_mod as tm

# Input file paths
turbie_parameters_file = "inputs/turbie_inputs/turbie_parameters.txt"
wind_file = "inputs/wind_files/wind_TI_0.1/wind_6_ms_TI_0.1.txt"
C_T_table_file = "inputs/turbie_inputs/CT.txt"

# Read the wind time series
time, wind_speed = tm.read_wind_file(wind_file)

# Read the thrust coefficient table
C_T_table = tm.read_C_T_table(C_T_table_file)

# Get thrust coefficient C_T from mean wind speed and C_T table
mean_wind_speed = np.mean(wind_speed)
C_T = tm.get_C_T(mean_wind_speed, C_T_table)

# Read Turbie parameters
parameters = tm.read_turbie_parameters(turbie_parameters_file)

# Get Turbie parameters
mb, mn, mh, mt, c1, c2, k1, k2, fb, ft, drb, drt, Dr, rho = tm.get_turbie_parameters(parameters)

# calculate M, C and K matrices
M, C, K = tm.calculate_MCK_matrix(mb, mn, mh, mt, c1, c2, k1, k2)

# Rotor area
A = tm.calculate_rotor_area(Dr)

# Test for the values of mean_wind_speed, C_T and A
print(mean_wind_speed, C_T, A)

# Initial state [x1, x2, x1_dot, x2_dot]: turbie start at rest
y0 = [0, 0, 0, 0]

# Solve the dynamics of the system using solve_ivp
sol = solve_ivp(
    tm.state_derivative,
    t_span=[time[0], time[-1]],
    y0=y0,
    t_eval=time,
    args=(M, C, K, time, wind_speed, C_T, rho, A)
)

# Check if the solution was successful
if sol.success:
    print("Simulation successful.", sol.y.shape)
    print("Mean blade displacement [m]:", np.mean(sol.y[0]))
    print("Mean tower displacement [m]:", np.mean(sol.y[1]))
else:
    print("Simulation failed.")

# Plot wind speed and displacements for this case
fig, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
fig.suptitle("Turbie response, mean wind speed 6 m/s, TI = 0.1")

# Panel 1: Wind speed vs time
axs[0].plot(time, wind_speed)
axs[0].set_ylabel("Wind Speed [m/s]")
axs[0].set_title("Wind Speed vs Time")

# Panel 2: Blade displacement x1 vs time
axs[1].plot(time, sol.y[0])
axs[1].set_ylabel("Blade Displacement [m]")
axs[1].set_title("Blade Displacement vs Time")

# Panel 3: Tower displacement x2 vs time
axs[2].plot(time, sol.y[1])
axs[2].set_ylabel("Tower Displacement [m]")
axs[2].set_title("Tower Displacement vs Time")
axs[2].set_xlim(time[0], time[-1])
axs[2].set_xlabel("Time [s]")

for ax in axs:
    ax.grid(True)
fig.tight_layout()

# Save the figure and show it
fig.savefig("outputs/time_series_6_ms_TI_0.1.png", dpi=300)
plt.show()