#!python3

"""
    Defining functions to calculate RMSE
"""

import math

import numpy as np
from scipy.interpolate import pchip_interpolate

"""
    Calculate population size of tennessen et al demography
        based on input time in the past
"""


def tennessen_popsize(t):
    if t < 204:
        return 5.12e6 * math.exp(-0.0195 * t)
    if t == 204:
        return 1.99e4
    if t > 204 and t < 920:
        return 1.99e4 * math.exp(-0.00307 * (t - 204))
    if t == 920:
        return 3.72e3
    if t > 920 and t < 2040:
        return 3.72e3
    if t == 2040:
        return 2.89e4
    if t > 2040 and t < 5920:
        return 2.89e4
    if t == 5920:
        return 1.46e4
    if t > 5920:
        return 1.46e4


"""
    Calculate RMSE
"""


def calc_rmse(y, y_truth):
    return np.sqrt(np.mean((y_truth - y) ** 2))


"""
    Calculate Mean Absolute Error
"""


def calc_mae(y, y_truth):
    return np.mean(np.abs(y_truth - y))


"""
    Calculates RMSE using spline interpolated values
"""


def calc_rmse_interpolated(t, y, lim, n):
    t_interp = np.linspace(0, lim, n)
    y_est_interp = pchip_interpolate(t, y, t_interp)
    Nt_true = np.array([tennessen_popsize(x) for x in t_interp])
    return calc_rmse(y_est_interp, Nt_true)
