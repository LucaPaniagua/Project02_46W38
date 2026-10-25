"""File turbie_mod.py"""
import numpy as np


# Read wind file and return the time and wind speed

def read_wind_file(filename):
    """Read a wind file and return the time and wind speed.

    Parameters
    ----------
    filename : str
        The path to the wind file.

    Returns
    -------
    time : np.ndarray
        An array containing the time values [s].
    wind_speed : np.ndarray
        An array containing the wind speed values [m/s].
    """

    time, wind_speed = np.loadtxt(filename, skiprows=1, unpack=True)
    return time, wind_speed

# C_T results from a mean wind speed and interpolation of the C_T table

def read_C_T_table(filename):
    """Read a C_T table from a file.

    Parameters
    ----------
    filename : str
        The path to the C_T table file.

    Returns
    -------
    C_T_table : np.ndarray
        A 2D array containing the wind speeds [m/s] and corresponding C_T values [-].
    """
    
    C_T_table = np.loadtxt(filename, skiprows=1)

    return C_T_table


def get_C_T(mean_wind_speed, C_T_table):
    """Get the thrust coefficient C_T from a mean wind speed and a C_T table.

    Parameters
    ----------
    mean_wind_speed : float
        The mean wind speed [m/s]
    C_T_table : np.ndarray
        A 2D array containing the wind speeds [m/s] and corresponding C_T values [-].

    Returns
    -------
    C_T : float
        The interpolated thrust coefficient [-].
    """
    
    wind_speeds = C_T_table[:, 0]
    C_T_values = C_T_table[:, 1]
    
    # Interpolate to find the corresponding C_T value for the given mean wind speed
    C_T = np.interp(mean_wind_speed, wind_speeds, C_T_values)
    
    return C_T



# Make matrix M, C  and K with turbie parameters

# Calculate the state derivate y´(t) for solve_ivp



