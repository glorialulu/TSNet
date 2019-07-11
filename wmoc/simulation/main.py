"""
The wmoc.simulation.main module contains function to perform 
the workflow of read, discretize, initial, and transient 
simulation for the given .inp file. 

"""
from __future__ import print_function
from wmoc.network import geometry, topology, discretization
from wmoc.network import valvesetting, pumpsetting
from wmoc.simulation.initialize import initialize
from wmoc.simulation.core import MOC
from datetime import datetime 
import numpy as np

def MOCSimulator(inp_file, dt, tf, valve_to_close=None, valve_op=None,
    pump_to_operate=None, pump_op=None,  
    leak_loc=None, leak_A=None, 
    burst_loc=None, burst_A=None, burst_t=None):
    r""" MOC Main Function
    
    Parameters
    ----------
    inp_file : .inp file
        Define the topology and parameters of the network   
    dt : float 
        Time step in transient simulation [s]
    tf : float 
        Duration of the transient simulation [s]
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
    ------
        wn : wntr.network.model.WaterNetworkModel
            Simulated network 
        H[npipe][nnode,tn] : numpy.array [m]
            Head results 
        V[npipe][nnode,tn] : numpy.array [m/s]
            Velocity results        
        tt : list [s]
            Simulated timestamps
    """

    startttime = datetime.now()
    # read the geometry of the network 
    wn, npipe = geometry(inp_file)
    # adjust the time step and discretize each pipe 
    wn, dt, Ndis = discretization(wn, npipe, dt)
    print('time step %.5f s' % dt)  
    
    tn = int(tf/dt) # Total time steps 
    print ('Total Time Step in this simulation', tn)

    # valve closure curve
    vo = valvesetting(dt, tf, valve_op)

    # pump operation curve
    po = pumpsetting(dt, tf, pump_op)

    # initial condition calculated at t0
    t0 = 0
    
    # determine leak location based on input node name
    # and add the leak to initial consition calculation
    if leak_loc != None :         
        for node in leak_loc:
            leak_node = wn.get_node(node)  
            leak_node.add_leak(wn, area=leak_A/np.sqrt(2*9.8)/1000, 
                                discharge_coeff=1,
                               start_time=t0)
        leak_loc = [wn.nodes[i].id-1 for i in leak_loc]

    # determine burst location based on input node name
    if burst_loc != None :
        burst_loc = [wn.nodes[i].id-1 for i in burst_loc]

    # calculate initial condition     
    wn, H, V = initialize(wn, t0, Ndis, npipe, tn)  

    # determine network topology 
    links1, links2, utype, dtype = topology(wn,npipe)    

    # startMOC transient simulation 
    H, V, tt = MOC(links1, links2, utype, dtype, wn,
                    H, V, npipe, Ndis, dt, tf, t0, 
                    valve_to_close, vo, pump_to_operate, po, 
                    leak_loc, leak_A, burst_loc, burst_A, burst_t)
    
    simtime = datetime.now() - startttime
    print('Running Time:', simtime)    
    return  wn, H, V, tt    
    