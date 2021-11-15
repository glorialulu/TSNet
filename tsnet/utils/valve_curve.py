"""
The tsnet.utils.valve_curve contains function to define
valve characteristics curve, gate valve by default.

"""
import numpy as np


def valve_curve(s, coeff=None):
    """Define valve curve

    Parameters
    ----------
    s : float
        open percentage
    valve : str, optional
        [description], by default 'Gate'

    Returns
    -------
    k : float
        Friction coefficient with given open percentage
    """
    if coeff is None:
        percent_open = np.linspace(100, 0, 11)
        # loss coefficients for a gate valve
        kl = [
            1 / 0.2, 2.50, 1.25, 0.625, 0.333, 0.17,
            0.100, 0.0556, 0.0313, 0.0167, 0.0
        ]
        k = np.interp(s, percent_open[::-1], kl[::-1])
    else:
        percent_open, kl = coeff
        k = np.interp(s, percent_open[::-1], kl[::-1])
    return k
