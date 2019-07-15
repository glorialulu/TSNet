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

def initialize(wn, t0, Ndis, npipe, tn):
    """Initial Condition Calculation.

    Initialize the list containing numpy arrays for velocity and head. 
    Calculate initial conditions using Epanet engine.
    Calculate D-W coefficients based on initial conditions. 
    Calculate demand coefficients based on initial conditions. 

    Parameters
    ----------
    wn : wntr.network.model.WaterNetworkModel
        Simulated network 
    t0 : float
        time to calculate initial condition
    npipe : float 
        Number of pipes 
    tn : float 
        Total time steps

    Returns
    -------
    wn : wntr.network.model.WaterNetworkModel
        Network with updated parameters 
    H : list 
        Initialized list of numpy.ndarray to store head results
    V : list
        Initialized list of numpy.ndarray to store velocity results
    """

    # initialize H and V for all computational nodes
    H = [0] * npipe  #pipes 
    V = [0] * npipe   
    
    # initialize the results list and set initial conditions   
    for _, pipe in wn.pipes():
        pn = int(pipe.id)-1
        # initialize the matrix for heads 
        H[pn] = np.zeros((int(Ndis[pn]+1),tn), dtype=np.float64)   
        # initialize the matrix for flow velocity                      
        V[pn] = np.zeros((int(Ndis[pn]+1),tn), dtype=np.float64) 

    print ("Initial Condition")            
    # calculate initial conditions using EPAnet engine
    sim = wntr.sim.EpanetSimulator(wn)
    results = sim.run_sim() 

    # assign the initial conditions to the result arrays
    for _, pipe in wn.pipes():      
        pn = int(pipe.id) -1 
        V[pn][:,0] = np.sign(results.link['flowrate'].loc[t0*3600, pipe.name])*\
                    results.link['velocity'].loc[t0, pipe.name]*\
                    np.ones((int(Ndis[pn]+1)))
        
        H[pn][:,0] = [results.node['head'].loc[t0, pipe.start_node_name] +\
                      i* ((results.node['head'].loc[t0, pipe.end_node_name]-
                       results.node['head'].loc[t0, pipe.start_node_name])/
                       int(Ndis[pn]+1)) 
                    for i in range(int(Ndis[pn]+1))]

        # calculate demand coefficient        
        Hs = H[pn][0,0]
        He = H[pn][-1,0]
        demand = [0,0]
        try :
            demand[0] = wn.nodes[pipe.start_node_name].demand_timeseries_list.at(t0)
        except:
            demand[0] = 0.
        try :
            demand[1] = wn.nodes[pipe.end_node_name].demand_timeseries_list.at(t0)
        except :
            demand[1] = 0.

        pipe = cal_demand_coef(demand, pipe, Hs, He, t0)

        # calculate roughness coefficien 
        Vp = V[pn][0,0]
        hl = abs(Hs - He)
        pipe = cal_roughness_coef(pipe, Vp, hl )
                 
    return wn, H, V

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

    pipe.start_demand_coeff = lambda: None # [m^3/s/(m H20)^(1/2)]
    pipe.end_demand_coeff = lambda: None   # [m^3/s/(m H20)^(1/2)]
    try:
        start_demand_coeff = demand[0]/ np.sqrt(Hs)
    except :
        start_demand_coeff = 0. 
        
    try: 
        end_demand_coeff = demand[1] / np.sqrt(He)
    except :
        end_demand_coeff = 0. 
        
    setattr(pipe, 'start_demand_coeff',start_demand_coeff )
    setattr(pipe, 'end_demand_coeff',end_demand_coeff )
    
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