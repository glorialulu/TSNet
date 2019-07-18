"""
The wmoc.simulation.main module contains function to perform
the workflow of read, discretize, initial, and transient
simulation for the given .inp file.

"""
from __future__ import print_function
from wmoc.network import  topology, discretization
from wmoc.network import valvesetting, pumpsetting
from wmoc.simulation.initialize import Initializer
from wmoc.simulation.core import MOC
from datetime import datetime
import numpy as np

def MOCSimulator(tm, valve_to_close=[], valve_op=None,
    pump_to_operate=[], pump_op=None,
    burst_loc=None, burst_A=None, burst_t=None):
    r""" MOC Main Function

    Parameters
    ----------
    tm : wmoc.network.model.TransientModel
        Network
    valve_to_close : list, optional
        The list of valves to be closed to generate transient, by default None
    valve_op : list, optional
        Required if valve_to_close is defined, by default None
        valve_op = [tc,ts,se,m]
        tc : the duration takes to close the valve [s]
        ts : closure start time [s]
        se : final open percentage [s]
        m  : closure constant [unitless]
    pump_to_operate: list, optional
        The list of pumps to be operate to generate transient, by default None
    pump_op : list, optional
        Required if pump_to_close is defined, by default None
        pump_op = [tc,ts,se,m]
        tc : the duration takes to operate the pump [s]
        ts : closure start time [s]
        se : final open percentage [s]
        m  : closure constant [unitless]
    burst_loc : list, optional
        The list of burst loction,
        defined by the name of the junction node, by default None
    burst_A : float, optional
        Required if burst_loc is defined
        The final burst coefficient of the burst, by default None
        Q_burst_final = burst_A [ m^3/s/(m H20)^(1/2)] * \sqrt(H)
    burst_t : list, optional
        Required if burst_loc is defined, by default None
        burst_t = [(burst_start_time, burst_end_time)] [s]

    Returns
    ------
    tm : wmoc.network.model.TransientModel
            Simulated network
    H[npipe][nnode,tn] : numpy.array [m]
        Head results
    V[npipe][nnode,tn] : numpy.array [m/s]
        Velocity results
    tt : list [s]
        Simulated timestamps
    """

    startttime = datetime.now()
    # initial condition calculated at t0
    t0 = 0

    # determine burst location based on input node name
    burst_A = None
    if burst_loc != None :
        burst_loc = [tm.nodes[i].id-1 for i in burst_loc]

    # determine network topology
    links1, links2, utype, dtype = topology(tm)

    # startMOC transient simulation
    tm = MOC(links1, links2, utype, dtype, tm, t0)

    simtime = datetime.now() - startttime
    print('Running Time:', simtime)
    return  tm
