import numpy as np
from scipy.interpolate import interp1d
import pandas as pd
from CoolProp.CoolProp import PropsSI
from math import pi, log


"""
Unit conversion functions
"""


def c2k(t):
    return t + 273.15


def lpm2cmps(v):
    return v * 1.6667e-5


"""
Get thermodynamic property functions
"""


def get_air_property_at1bar(which_property, which_temp):
    """
    This function is interpolating some thermodynamic properties of air at 1 bar.
    Boundary values taken from one of Cengel's books.
    Created for manual calculation of the surface heat transfer coefficient (radiation + forced + nat convection)
    :param which_property: property name to be interpolated
    :param which_temp: temperature of the desired property
    :return: property value interpolated for the desired temperature
    """
    air_props_header = ['rho', 'cp', 'k', 'alpha', 'mu', 'nu', 'pr']
    t10 = [1.246, 1006, 0.02439, 1.944e-5, 1.778e-5, 1.426e-5, 0.7336]
    # t20 = [1.204, 1007, 0.02514, 2.074e-5, 1.825e-5, 1.516e-5, 0.7309]
    t30 = [1.164, 1007, 0.02588, 2.208e-5, 1.872e-5, 1.608e-5, 0.7282]
    zipped = zip(t10, t30)
    df = pd.DataFrame(zipped, index=air_props_header, columns=['deg10', 'deg30'])
    pd.set_option('display.float_format', '{:.5g}'.format)
    prop_interp_func = interp1d([10, 30], [df.deg10[which_property], df.deg30[which_property]])
    new_prop = prop_interp_func(which_temp)
    return new_prop


def get_cp(t, fluid='IF97::Water', p=101325):
    return PropsSI('C', 'T', t, 'P', p, fluid)


def get_rho(t, fluid='IF97::Water', p=101325):
    return PropsSI('D', 'T', t, 'P', p, fluid)


def get_k(material, t):
    if material is 'tubolit':
        return (36 + 0.1 * t + 0.0008 * (t-40)**2) / 1000
    elif material is 'armaflexaf':
        return (33 + 0.1 * t + 0.0008 * t ** 2) / 1000
    elif material is 'copper':
        return 397
    else:
        return print('No such material defined yet')


"""
Other functions
"""


def get_rmse(pred, obs):
    return int(round(np.sqrt(((pred - obs)**2).mean())))


def get_t_avg(t_in, t_out):
    return (t_in + t_out) / 2


def get_surf_area_cyl(radius, length):
    return 2 * pi * radius * length


def get_cond_res(r_outer, r_inner, conductivity, length):
    return log(r_outer / r_inner) / (2 * pi * conductivity * length)
