"""
The wmoc.simulation.initialize contains functions to
1. Initialize the list containing numpy arrays for velocity and head.
2. Calculate initial conditions using Epanet engine.
3. Calculate D-W coefficients based on initial conditions.
4. Calculate demand coefficients based on initial conditions.

"""

import wntr
import numpy as np
import warnings
from wmoc.network import discretization

def Initializer(tm, t0, dt, tf):
    """Initial Condition Calculation.

    Initialize the list containing numpy arrays for velocity and head.
    Calculate initial conditions using Epanet engine.
    Calculate D-W coefficients based on initial conditions.
    Calculate demand coefficients based on initial conditions.

    Parameters
    ----------
    tm : wmoc.network.geometry.TransientModel
        Simulated network
    t0 : float
        time to calculate initial condition
    dt : float
        User defined time step
    tf : float or int
        Simulation duration

    Returns
    -------
    tm : wmoc.network.geometry.TransientModel
        Network with updated parameters
    H : list
        Initialized list of numpy.ndarray to store latest head results
    V : list
        Initialized list of numpy.ndarray to store latest velocity results
    """
    # adjust the time step and discretize each pipe
    tm = discretization(tm, dt)
    print('Simulation time step %.5f s' % tm.time_step)
    tm.simulation_peroid = tf
    tn = int(tf/tm.time_step) # Total time steps
    print ('Total Time Step in this simulation %s'  %tn)

    # create new atributes for each pipe to store head and velocity results
    # at its start and end node.

    for _, pipe in tm.pipes():
        pipe.start_node_head = np.zeros(tn)
        pipe.start_node_velocity = np.zeros(tn)
        pipe.end_node_head = np.zeros(tn)
        pipe.end_node_velocity = np.zeros(tn)

    # calculate initial conditions using EPAnet engine
    sim = wntr.sim.WNTRSimulator(tm)
    results = sim.run_sim()

    for _, pipe in tm.pipes():
        # assign the initial conditions to the latest result arrays

        V = np.sign(results.link['flowrate'].loc[t0, pipe.name])*\
                    results.link['velocity'].loc[t0, pipe.name]*\
                    np.ones(pipe.number_of_segments+1)

        H = [results.node['head'].loc[t0, pipe.start_node_name] +\
                      i* ((results.node['head'].loc[t0, pipe.end_node_name]-
                       results.node['head'].loc[t0, pipe.start_node_name])/
                       pipe.number_of_segments)
                    for i in range(pipe.number_of_segments+1)]
                    
        H = np.array(H)
        pipe.initial_head = H
        pipe.initial_velocity = V

        # assign the initial conditions to the results attributes
        pipe.start_node_velocity[0] = V[0]
        pipe.end_node_velocity[0] = V[-1]
        pipe.start_node_head[0] = H[0]
        pipe.end_node_head[0] = H[-1]
        # print('test', pipe.name, pipe.start_node_velocity[0],V[0])
        # calculate demand coefficient
        Hs = H[0]
        He = H[-1]
        demand = [0,0] # dmenad at start and end node

        try :
            demand[0] = tm.nodes[pipe.start_node_name].demand_timeseries_list.at(t0)
        except:
            demand[0] = 0.
        try :
            demand[1] = tm.nodes[pipe.end_node_name].demand_timeseries_list.at(t0)
        except :
            demand[1] = 0.

        pipe = cal_demand_coef(demand, pipe, Hs, He, t0)

        # calculate roughness coefficient
        Vp = V[0]
        hl = abs(Hs - He)
        pipe = cal_roughness_coef(pipe, Vp, hl )

    # set initial conditions as a new attribute to TransientModel
    tm.initial_head = H
    tm.initial_velocity = V

    return tm

def cal_demand_coef(demand, pipe, Hs, He, t0=0.):

    """Calculate the demand coefficent for the start and end node of the pipe.

    Parameters
    ----------
    demand : list
        Demand at the start (demand[0]) and end demand[1] node
    pipe : object
        Pipe object
    Hs : float
        Head at the start node
    He : float
        Head at the end node
    t0 : float, optional
        Time to start initial condition calculation, by default 0

    Returns
    -------
    pipe : object
        Pipe object with calculated demand coefficient
    """


    try:
        start_demand_coeff = demand[0]/ np.sqrt(Hs)
    except :
        start_demand_coeff = 0.

    try:
        end_demand_coeff = demand[1] / np.sqrt(He)
    except :
        end_demand_coeff = 0.
    pipe.start_demand_coeff = start_demand_coeff # [m^3/s/(m H20)^(1/2)]
    pipe.end_demand_coeff = end_demand_coeff   # [m^3/s/(m H20)^(1/2)]


    return pipe

def cal_roughness_coef(pipe, V, hl):
    """Calculate the D-W roughness coefficient based on initial conditions.

    Parameters
    ----------
    pipe : object
        Pipe object
    V : float
        Initial flow velocity in the pipe
    hl : float
        Initial head loss in the pipe

    Returns
    -------
    pipe : object
        Pipe object with calculated D-W roughness coefficient.
    """
    g = 9.8
    H_tol = 1e-4
    V_tol = 1e-5

    if abs(V) >= V_tol and hl >= H_tol:
        pipe.roughness = hl / (pipe.length/pipe.diameter) / (V**2/2/g)
    else:
        pipe.roughness = 0

    if pipe.roughness >0.08:
        warnings.warn("%s :the friction coefficient %.4f is too large. \
                        The D-W coeff has been set to 0.03 "
                        %(pipe.name, pipe.roughness))
        pipe.roughness = 0.03

    return pipe