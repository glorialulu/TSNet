"""
The tsnet.postprocessing package contains functions to
postprocess the simulation results.

"""
from .time_history import plot_head_history, plot_velocity_history
from .detect_cusum import detect_cusum

__all__ = [
    'plot_head_history',
    'plot_velocity_history',
    'detect_cusum'
]
