"""
The wmoc.network package contains methods to define 
1. a water network geometry,
2. network topology, 
3. network control, and
4 .spatial and temporal discretization.

"""

from .topology import topology
from .geometry import TransientModel
from .discretize import discretization
from .control import valvesetting, pumpsetting
