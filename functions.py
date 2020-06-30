import numpy as np

"""
Unit conversion functions
"""


def c2k(t):
    return t + 273.15


def lpm2cmps(v):
    return v * 1.6667e-5


"""
Other functions
"""


def get_rmse(pred, obs):
    return int(round(np.sqrt(((pred - obs)**2).mean())))