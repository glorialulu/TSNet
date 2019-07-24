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

def Initializer(tm, t0, engine='PDD'):
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
    engine : string
        steady state calculation engine:
        DD: demand driven;
        PDD: pressure dependent demand,
        by default DD


    Returns
    -------
    tm : wmoc.network.geometry.TransientModel
        Network with updated parameters
    """
    # adjust the time step and discretize each pipe

    tn = int(tm.simulation_peroid/tm.time_step) # Total time steps
    print ('Total Time Step in this simulation %s'  %tn)

    # create new attributes for each pipe to store head and velocity results
    # at its start and end node.
    for _, pipe in tm.pipes():
        pipe.start_node_head = np.zeros(tn)
        pipe.start_node_velocity = np.zeros(tn)
        pipe.start_node_flowrate = np.zeros(tn)
        pipe.end_node_head = np.zeros(tn)
        pipe.end_node_velocity = np.zeros(tn)
        pipe.end_node_flowrate = np.zeros(tn)

    # create new atributes for each node to store head and discharge results
    for _,node in tm.nodes():
        node.demand_discharge = np.zeros(tn)
        node.emitter_discharge = np.zeros(tn)

    # calculate initial conditions using EPAnet engine
    for _,node in tm.nodes():
        if node.leak_status == True:
            node.add_leak(tm, area=node.emitter_coeff/np.sqrt(2*9.8),
                    discharge_coeff = 1, start_time = t0)
    if engine.lower() == 'dd':
        sim = wntr.sim.WNTRSimulator(tm)
        results = sim.run_sim()
    elif engine.lower() == 'pdd':
        sim = wntr.sim.WNTRSimulator(tm,mode='PDD')
        results = sim.run_sim()
    else:
        raise Exception("Unknown initial calculation engine. \
            The engine can only be 'WNTR' or 'EPANET'.")

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
        pipe.start_node_flowrate[0] = V[0]*pipe.area
        pipe.end_node_flowrate[0] = V[-1]*pipe.area

        try:
            pipe.start_node.emitter_discharge[0] = pipe.start_node.emitter_coeff * np.sqrt(H[0])
            pipe.start_node.demand_discharge[0] = pipe.start_node.demand_coeff * np.sqrt(H[0])
        except:
            pass # reservoir does not have emitter_discharge attribute

        try:
            pipe.end_node.emitter_discharge[0] = pipe.end_node.emitter_coeff * np.sqrt(H[-1])
            pipe.end_node.demand_discharge[0] =  pipe.end_node.demand_coeff * np.sqrt(H[-1])
        except:
            pass # reservoir does not have emitter_discharge attribute

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
    pipe.start_node.demand_coeff = start_demand_coeff # [m^3/s/(m H20)^(1/2)]
    pipe.end_node.demand_coeff = end_demand_coeff   # [m^3/s/(m H20)^(1/2)]


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
    H_tol = 1e-3
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