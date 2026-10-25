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

def read_turbie_parameters(filename):
    """Read the Turbie parameters from a file and return an array with the parameters.

    Parameters
    ----------
    filename : str
        The path to the Turbie parameters file.

    Returns
    -------
    parameters : np.ndarray
        1D array with the 14 Turbie parameters, in file order:
        mb, mn, mh, mt [kg], c1, c2 [N/(m/s)], k1, k2 [N/m],
        fb, ft [Hz], drb, drt [-], Dr [m], rho [kg/m^3].
    """

    parameters = np.loadtxt(filename, skiprows=1)

    return parameters

def get_turbie_parameters(parameters):
    """Extract the Turbie parameters from the array and return them as individual variables.

    Parameters
    ----------
    parameters : np.ndarray
        1D array with the 14 Turbie parameters, in file order:
        mb, mn, mh, mt [kg], c1, c2 [N/(m/s)], k1, k2 [N/m],
        fb, ft [Hz], drb, drt [-], Dr [m], rho [kg/m^3].

    Returns
    -------
    mb : float
        Mass of the base [kg].
    mn : float
        Mass of the nacelle [kg].
    mh : float
        Mass of the hub [kg].
    mt : float
        Mass of the tower [kg].
    c1 : float
        Damping coefficient for the base [N/(m/s)].
    c2 : float
        Damping coefficient for the tower [N/(m/s)].
    k1 : float
        Stiffness coefficient for the base [N/m].
    k2 : float
        Stiffness coefficient for the tower [N/m].
    fb : float
        Natural frequency of the base [Hz].
    ft : float
        Natural frequency of the tower [Hz].
    drb : float
        Damping ratio for the base [-].
    drt : float
        Damping ratio for the tower [-].
    Dr : float
        Rotor diameter [m].
    rho : float
        Air density [kg/m^3].
    """
    
    mb, mn, mh, mt = parameters[0:4]
    c1, c2 = parameters[4:6]
    k1, k2 = parameters[6:8]
    fb, ft = parameters[8:10]
    drb, drt = parameters[10:12]
    Dr = parameters[12]
    rho = parameters[13]

    return mb, mn, mh, mt, c1, c2, k1, k2, fb, ft, drb, drt, Dr, rho

def calculate_MCK_matrix(mb, mn, mh, mt, c1, c2, k1, k2):
    """Calculate the mass, damping, and stiffness matrices for the Turbie system.

    Parameters
    ----------
    mb : float
        Mass of the base [kg].
    mn : float
        Mass of the nacelle [kg].
    mh : float
        Mass of the hub [kg].
    mt : float
        Mass of the tower [kg].
    c1 : float
        Damping coefficient for the base [N/(m/s)].
    c2 : float
        Damping coefficient for the tower [N/(m/s)].
    k1 : float
        Stiffness coefficient for the base [N/m].
    k2 : float
        Stiffness coefficient for the tower [N/m].

    Returns
    -------
    M : np.ndarray
        The mass matrix. m1 = 3*mb; m2 = mn + mh + mt
    C : np.ndarray
        The damping matrix.
    K : np.ndarray
        The stiffness matrix.
    """

    m1 = 3*mb
    m2 = mn + mh + mt
    
    M = np.array([[m1, 0], [0, m2]])
    C = np.array([[c1, -c1], [-c1, c1+c2]])
    K = np.array([[k1, -k1], [-k1, k1+k2]])
    
    return M, C, K


# Calculate the state derivate y´(t) for solve_ivp





