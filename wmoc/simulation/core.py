"""
The wmoc.simulation.core module contains the core function to
perform trasnient simulation.

"""
from __future__ import print_function
from wmoc.simulation.single import inner_pipe, left_boundary, right_boundary
from wmoc.utils import valve_curve
import numpy as np
from functools import update_wrapper
import warnings

def decorator(d):
    """Make function d a decorator: d wraps a function fn."""
    def _d(fn):
        return update_wrapper(d(fn), fn)
    update_wrapper(_d, d)
    return _d

def memo(f):
    """Decorator that caches the return value for each call to f(args).
    Then when called again with same args, we can just look it up."""
    cache = {}
    def _f(*args):
        try:
            return cache[args]
        except KeyError:
            cache[args] = result = f(*args)
            return result
        except TypeError:
            # some element of args can't be a dict key
            return f(*args)
    return _f


def MOC(links1, links2, utype, dtype, tm, t0,
    burst_loc=None, final_burst=None, burst_t=None):
    r"""Transient Simulation using MOC method

    The core function to perform transient simulation
    using MOC method.

    Parameters
    ----------
    links1 : list
        The id of ajacent pipe on the start node.
        The sign represents the direction of the pipe.
        + : flowing into the junction
        - : flowing out from the junction
    links2 : list
        The id of ajacent pipe on the end node.
        The sign represents the direction of the pipe.
        + : flowing into the junction
        - : flowing out from the junction
    utype : list
        The type of the upstream ajacent links.
        If the link is not pipe, the name of that link
        will also be included.
        If there is no upstream link, the type of the start node
        will be recorded.
    dtype : list
        The type of the downstream ajacent links.
        If the link is not pipe, the name of that link
        will also be included.
        If there is no downstream link, the type of the end node
        will be recorded.
    tm : wmoc.network.geometry.TransientModel
        Network
    H : list
        Array to store head results with initial condition given
    V : list
        Array to store head results with initial condition given
    npipe : integer
        Number of pipes
    Ndis : numpy array
        Number of discritization for each pipe
    dt : float
        Adjusted time step
    tf : float
        Simulation Time
    t0 : float
        Time to calculate initial condition
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
    leak_loc : list, optional
        The list of leakage loction, defined by the name of
        the junction node, by default None
    leak_A : float, optional
        Required if leak_loc is defined
        The leakage coefficient of the leakge, by default None
        .. math:: Q_leak = leak_A  [ m^3/s/(m H20)^(1/2)] * \sqrt(H)
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
    -------
    H[npipe][nnode,tn] : numpy.ndarray [m]
        Head results
    V[npipe][nnode,tn] : numpy.ndarray [m/s]
        Velocity results
    tt : list [s]
        Simulated timestamps
    """

    leak_loc = [name.id-1 for name in tm.leak_node]
    leak_A = [node.emitter_coeff for node in tm.leak_node]
    tt = ['x']
    tt.append(0)
    dt = tm.time_step
    tn = int(tm.simulation_peroid/tm.time_step)  # Total time steps

    # determine which node of the adjacant pipe should be call:
    # if the adjacant pipe is entering the junction, then -2
    # if the adjacent pipe is leaving the junction, then 1
    a = {1:-2, -1:1}
    # generat a list of pipe
    p = []
    # results from last time step
    H = [0] * tm.num_pipes
    V = [0] * tm.num_pipes
    # results at current time step
    HN = [0] * tm.num_pipes
    VN = [0] * tm.num_pipes

    for _, pipe in tm.pipes():
        p.append(pipe)

    # initial condition
    for _, pipe in tm.pipes():
        pn = pipe.id-1
        H[pn] = pipe.initial_head
        V[pn] = pipe.initial_velocity

    # Start Claculation
    for ts in range(1,tn):
        t = ts*dt
        tt.append(t)
        # initilaize the results at this time step
        for _, pipe in tm.pipes():
            pn = pipe.id-1
            HN[pn] =  np.zeros_like(H[pn])
            VN[pn] =  np.zeros_like(V[pn])

        # calculate the burst area in the current time
        burst_A = None
        if burst_loc != None:
            if (burst_t[0][1]-burst_t[0][0]) !=0:
                tmp = (t - burst_t[0][0])/(burst_t[0][1]-burst_t[0][0])
                if tmp < 0:
                    burst_A = 0
                elif tmp > 1 :
                    burst_A = final_burst
                else :
                    burst_A = final_burst * tmp
            else:
                if t < burst_t[0][0]:
                    burst_A = 0
                else :
                    burst_A = final_burst

        for _,pipe in tm.pipes():
            pn = pipe.id-1
            # Assumption:
            # when a pipe is connected with a pump or valve,
            # the connection is not branch junction.

            # inner pipes
            if links1[pn] and links2[pn] and \
                links1[pn] != ['End'] and links2[pn] != ['End']:
                # list to store information about pump and vale
                # pump[0] and valve[0] for upstream elemnets
                # pump[1] and valve[1] for downstream elements
                pump = [[],[]]; valve = [0,0]
                # upstream
                if utype[pn][0] == 'Pump':
                    # three points for pump charatersitics curve
                    pump[0] = tm.links[utype[pn][1]].get_pump_curve().points
                    # calculate the coordinate of the three points
                    # based on the pump speed
                    if tm.links[utype[pn][1]].operating == True:
                        po = tm.links[utype[pn][1]].operation_rule[ts]
                        pump[0]=[(i*po,j*po**2) for (i,j) in pump[0]]

                elif utype[pn][0] == 'Valve':
                    # determine valve fricton coefficients based on
                    # open percentage
                    if tm.links[utype[pn][1]].operating == True:
                        valve[0] = valve_curve(tm.links[utype[pn][1]].operation_rule[ts]*100)
                    else :
                         valve[0] = valve_curve(100)
                # downstream
                if dtype[pn][0] == 'Pump':
                    pump[1] = tm.links[dtype[pn][1]].get_pump_curve().points
                    if tm.links[dtype[pn][1]].operating == True:
                        po = tm.links[dtype[pn][1]].operation_rule[ts]
                        pump[1]=[(i*po,j*po**2) for (i,j) in pump[1]]

                elif dtype[pn][0] == 'Valve':
                    if tm.links[dtype[pn][1]].operating == True:
                        valve[1] = valve_curve(tm.links[dtype[pn][0]].operation_rule[ts]*100)
                    else :
                         valve[1] = valve_curve(100)

                HN[pn], VN[pn] = inner_pipe(pipe, pn, dt,
                     links1[pn], links2[pn], utype[pn], dtype[pn], p,
                     H[pn], V[pn], HN[pn], VN[pn],
                     [H[abs(i)-1][a[np.sign(i)]] for i in links1[pn]],
                     [V[abs(i)-1][a[np.sign(i)]] for i in links1[pn]],
                     [H[abs(i)-1][a[np.sign(i)]] for i in links2[pn]],
                     [V[abs(i)-1][a[np.sign(i)]] for i in links2[pn]],
                     pump, valve,
                     leak_loc, leak_A, burst_loc, burst_A)

                pipe.start_node_velocity[ts] = VN[pn][0]
                pipe.end_node_velocity[ts] = VN[pn][-1]
                pipe.start_node_head[ts] = HN[pn][0]
                pipe.end_node_head[ts] = HN[pn][-1]

            # left boundary pipe
            elif not links1[pn] or links1[pn] == ['End']:
                pump = [[],[]]; valve = [0,0]
                # LEFT BOUNDARY
                if utype[pn][0] == 'Reservoir':
                    # head B.C.
                    HN[pn][0]   =  tm.nodes[utype[pn][1]].base_head
                elif utype[pn][0] == 'Tank':
                    # head B.C.
                    HN[pn][0]   =  tm.nodes[utype[pn][1]].head
                elif utype[pn][0] == 'Junction':
                    VN[pn][0] = pipe.initial_velocity[0]
                elif utype[pn][0] == 'Valve':
                    if tm.links[utype[pn][1]].operating == True:
                        # velocity B.C.
                        VN[pn][0]   = pipe.initial_velocity[0] * \
                            tm.links[utype[pn][1]].operation_rule[ts]
                    else :
                        VN[pn][0]   = pipe.initial_velocity[0]
                # source pump
                elif utype[pn][0] == 'Pump':
                    # pump[0][0]: elevation of the reservoir/tank
                    # pump[0][1]: three points for pump characteristic curve
                    pump[0] = [tm.links[utype[pn][1]].start_node.base_head,
                         tm.links[utype[pn][1]].get_pump_curve().points]

                    if tm.links[utype[pn][1]].operating == True:
                        po = tm.links[utype[pn][1]].operation_rule[ts]
                        pump[0][1]=[(i*po,j*po**2) for (i,j) in pump[0][1]]
                        print('Operating pump %s' %utype[pn][1])
                else:
                     warnings.warn ('Pipe %s miss %s upstream.' %(pipe, utype[pn][0]))

                # RIGHT BOUNDARY
                if dtype[pn][0] == 'Pump':
                    pump[1] = tm.links[dtype[pn][1]].get_pump_curve().points
                    if tm.links[dtype[pn][1]].operating == True:
                        po = tm.links[dtype[pn][1]].operation_rule[ts]
                        pump[1]=[(i*po[ts],j*po[ts]**2) for (i,j) in pump[1]]

                elif dtype[pn][0] == 'Valve':
                    if tm.links[dtype[pn][1]].operating == True:
                        valve[1] = valve_curve(tm.links[dtype[pn][1]].operation_rule*100)
                    else :
                         valve[1] = valve_curve(100)

                HN[pn], VN[pn] = left_boundary(pipe, pn,
                      HN[pn], VN[pn], H[pn], V[pn],
                     links2[pn], p, pump, valve, dt,
                     [H[abs(i)-1][a[np.sign(i)]] for i in links2[pn]],
                     [V[abs(i)-1][a[np.sign(i)]] for i in links2[pn]],
                     utype[pn], dtype[pn], leak_loc, leak_A, burst_loc, burst_A)

                pipe.start_node_velocity[ts] = VN[pn][0]
                pipe.end_node_velocity[ts] = VN[pn][-1]
                pipe.start_node_head[ts] = HN[pn][0]
                pipe.end_node_head[ts] = HN[pn][-1]

            #  right boundary pipe
            elif not links2[pn] or links2[pn] == ['End']:
                pump = [[],[]]; valve = [0,0]
                # RIGHT boundary
                if dtype[pn][0] == 'Reservoir':
                    HN[pn][-1]   =  tm.nodes[dtype[pn][1]].base_head # head of reservoir
                elif dtype[pn][0] == 'Tank':
                    HN[pn][-1]   =  tm.nodes[dtype[pn][1]].head # head of tank
                elif dtype[pn][0] == 'Junction':
                    VN[pn][-1] = pipe.initial_velocity[-1]
                elif dtype[pn][0] == 'Valve':
                    if tm.links[dtype[pn][1]].operating == True:
                        # valve velocity condition
                        VN[pn][-1] = pipe.initial_velocity[-1]* \
                        tm.links[dtype[pn][1]].operation_rule[ts]  
                    else :
                        VN[pn][-1] = pipe.initial_velocity[-1]
                else :
                     warnings.warn('Pipe %s miss %s downstream.' %(pipe, dtype[pn][0]))
                # LEFT boundary
                # source pump
                if utype[pn][0] == 'Pump':
                    # pump[0][0]: elevation of the reservoir/tank
                    # pump[0][1]: three points for pump characteristic curve
                    pump[0] = tm.links[utype[pn][1]].get_pump_curve().points
                    if tm.links[utype[pn][1]].operating == True:
                        po = tm.links[utype[pn][1]].operation_rule[ts]
                        pump[0]=[(i*po,j*po**2) for (i,j) in pump[0]]

                elif utype[pn][0] == 'Valve':
                    if tm.links[utype[pn][1]].operating == True:
                        valve[0] = valve_curve(tm.links[utype[pn][1]].operation_rule[ts]*100)
                    else :
                         valve[0] = valve_curve(100)

                HN[pn], VN[pn] = right_boundary(pipe, pn,
                     H[pn], V[pn], HN[pn], VN[pn],
                     links1[pn], p, pump, valve,  dt,
                     [H[abs(i)-1][a[np.sign(i)]] for i in links1[pn]],
                     [V[abs(i)-1][a[np.sign(i)]] for i in links1[pn]],
                     utype[pn], dtype[pn], leak_loc, leak_A, burst_loc, burst_A)

                pipe.start_node_velocity[ts] = VN[pn][0]
                pipe.end_node_velocity[ts] = VN[pn][-1]
                pipe.start_node_head[ts] = HN[pn][0]
                pipe.end_node_head[ts] = HN[pn][-1]

        # march in time
        for _, pipe in tm.pipes():
            pn = pipe.id-1
            H[pn] = HN[pn]
            V[pn] = VN[pn]

    tm.simulation_timestamps = tt[1:]
    return tm