"""
The wmoc.simulation.single contains methods to perform MOC 
transient simulation on a single pipe, including 
1. inner pipe
2. left boundary pipe (without C- charateristic grid)
3. right boundary pipr (without C+ characteristic grid)

"""
import numpy as np 
from WMOC.simulation.solver import *
from WMOC.utils import calc_parabola_vertex

def inner_pipe (linkp, pn,  H, V, links1, links2, utype, dtype, p, n, dt, 
                H0, V0, H10, V10, H20, V20, pump, valve,
                leak_loc, leak_A, burst_loc, burst_A):
    """MOC solution for an indivial inner pipe.
    
    Parameters
    ----------
    linkp : object
        Current pipe object
    pn : int
        Current pipe ID
    H : numpy.ndarray 
        Head of current pipe at current time step [m]
    V : numpy.ndarray
        Velocity of current pipe at current time step [m/s]
    links1 : list
        Upstream adjacent pipes 
    links2 : list
        Downstream ajacent pipes 
    utype : list 
        Upstream adjacent link type, and if not pipe, their name
    dtype : list
        Downstream adjacent link type, and if not pipe, their name
    p : list
        pipe list 
    n : int
        Number of discretization of current pipe
    dt : float 
        Time step 
    H0 : numpy.ndarray
        Head of current pipe at previous time step [m]
    V0 : numpy.ndarray
        Velocity of current pipe at previous time step [m/s]
    H10 : list
        Head of left adjacent nodes at previous time step [m]
    V10 : list
        Velocity of left adjacent nodes at previous time step [m/s]
    H20 : list
        Head of right adjacent nodes at previous time step [m]
    V20 : list
        Velocity of right adjacent nodes at previous time step [m/s]
    pump : list
        Characteristics of the pump 
    valve : list
        Characteristics of the valve
    leak_loc : list, optional 
        The list of leakage loction, defined by the name of
        the junction node, by default None
    leak_A : float, optional
        Required if leak_loc is defined
        The leakage coefficient of the leakge, by default None 
        Q_leak = leak_A  [ m^3/s/(m H20)^(1/2)] * \sqrt(H)
    burst_loc : list, optional 
        The list of burst loction,
        defined by the name of the junction node, by default None
    burst_A : float, optional 
        Required if burst_loc is defined
        The burst coefficient of the burst, by default None 
        Q_burst_final = burst_A [ m^3/s/(m H20)^(1/2)] * \sqrt(H)
    
    Returns
    -------
    H : numpy.ndarray
        Head results of the current pipe at current time step. [m]
    V : numpy.ndarray
        Velocity results of the current pipe at current time step. [m/s]
    """

    # Properties of current pipe 
    g = 9.8                          # m/s^2
    n = int(n)                       # spatial discritization 
    link1 = [p[abs(i)-1] for i in links1]
    link2 = [p[abs(i)-1] for i in links2]
    
    for i in range(n+1): 
        # Pipe start
        if i == 0:
            V1 = V10;     H1 = H10       #list
            V2 = V0[i+1]; H2 = H0[i+1]
            if utype[0] == 'Pipe':
                if leak_loc != None and int(p[pn].start_node.id) in leak_loc:
                    leak_A += linkp.start_demand_coeff
                    H[i], V[i] = add_leakage(leak_A, link1, linkp,
                     H1, V1, H2, V2, dt, g, i,  np.sign(links1), [-1])
                    
                elif burst_loc != None and int(p[pn].start_node.id) in burst_loc:
                    burst_A += linkp.start_demand_coeff
                    H[i], V[i] = add_leakage(burst_A, link1, linkp,
                     H1, V1, H2, V2, dt, g, i,  np.sign(links1), [-1])
                    
                else :
                    H[i], V[i] = add_leakage(linkp.start_demand_coeff, link1, linkp,
                     H1, V1, H2, V2, dt, g, i,  np.sign(links1), [-1])
            elif utype[0] == 'Pump':
                pumpc = calc_parabola_vertex(pump[0])
                H[i], V[i] = pump_node(pumpc, link1, linkp,
                     H1, V1, H2, V2, dt, g, i,  np.sign(links1), [-1])
            elif utype[0] == 'Valve':
                valvec = valve[0]
                H[i], V[i] = valve_node(valvec, link1, linkp,
                     H1, V1, H2, V2, dt, g, i,  np.sign(links1), [-1])
        # Pipe end  
        if i == n:
            V1 = V0[i-1];    H1 = H0[i-1]
            V2 = V20;        H2 = H20  
            if dtype[0] == 'Pipe':
                if leak_loc != None and int(p[pn].end_node.id) in leak_loc:
                    leak_A += linkp.end_demand_coeff
                    H[i], V[i] = add_leakage(leak_A,linkp, link2,
                     H1, V1, H2, V2, dt, g, i, [1], np.sign(links2))
                elif burst_loc != None and int(p[pn].end_node.id) in burst_loc:
                    burst_A += linkp.end_demand_coeff
                    H[i], V[i] = add_leakage(burst_A,linkp, link2,
                     H1, V1, H2, V2, dt, g, i, [1], np.sign(links2))
                else :
                    H[i], V[i] = add_leakage(linkp.end_demand_coeff,linkp, link2,
                     H1, V1, H2, V2, dt, g, i, [1], np.sign(links2))
                    

                    
            elif dtype[0] == 'Pump':
                pumpc = calc_parabola_vertex(pump[1])
                H[i], V[i] = pump_node(pumpc, linkp, link2,
                     H1, V1, H2, V2, dt, g, i, [1], np.sign(links2))
                
            elif dtype[0] == 'Valve':
                valvec = valve[1]
                H[i], V[i] = valve_node(valvec, linkp, link2,
                     H1, V1, H2, V2, dt, g, i, [1], np.sign(links2))
                
        # Interior points
        if (i > 0) and (i < n):
            V1 = V0[i-1]; H1 = H0[i-1]
            V2 = V0[i+1]; H2 = H0[i+1]
                
            H[i], V[i] = inner_node(linkp, linkp, 0,
             H1, V1, H2, V2, dt, g, i,[1],[-1])
                  
    return H, V 

def left_boundary(linkp, pn, H, V, links2, p, pump, valve, n, dt, 
                  H0, V0, H20, V20, utype, dtype,
                  leak_loc, leak_A, burst_loc, burst_A) :
    """MOC solution for an indivial left boundary pipe.
    
    Parameters
    ----------
    linkp : object
        Current pipe object
    pn : int
        Current pipe ID
    H : numpy.ndarray 
        Head of current pipe at current time step [m]
    V : numpy.ndarray
        Velocity of current pipe at current time step [m/s]
    links2 : list
        Downstream ajacent pipes 
    p : list
        pipe list 
    pump : list
        Characteristics of the pump 
    valve : list
        Characteristics of the valve
    n : int
        Number of discretization of current pipe
    dt : float 
        Time step 
    H0 : numpy.ndarray
        Head of current pipe at previous time step [m]
    V0 : numpy.ndarray
        Velocity of current pipe at previous time step [m/s]
    H20 : list
        Head of right adjacent nodes at previous time step [m]
    V20 : list
        Velocity of right adjacent nodes at previous time step [m/s]
    utype : list 
        Upstream adjacent link type, and if not pipe, their name
    dtype : list
        Downstream adjacent link type, and if not pipe, their name
    leak_loc : list, optional 
        The list of leakage loction, defined by the name of
        the junction node, by default None
    leak_A : float, optional
        Required if leak_loc is defined
        The leakage coefficient of the leakge, by default None 
        Q_leak = leak_A  [ m^3/s/(m H20)^(1/2)] * \sqrt(H)
    burst_loc : list, optional 
        The list of burst loction,
        defined by the name of the junction node, by default None
    burst_A : float, optional 
        Required if burst_loc is defined
        The burst coefficient of the burst, by default None 
        Q_burst_final = burst_A [ m^3/s/(m H20)^(1/2)] * \sqrt(H)
    
    Returns
    -------
    H : numpy.ndarray
        Head results of the current pipe at current time step. [m]
    V : numpy.ndarray
        Velocity results of the current pipe at current time step. [m/s]
    """
    
    link2 = [p[abs(i)-1] for i in links2]
    # Properties of current pipe 
    f = linkp.roughness              # unitless
    D = linkp.diameter               # m
    g = 9.8                          # m/s^2
    a = linkp.wavev    # m/s
    n = int(n)                       # spatial discritization 
      
    for i in range(n+1):
        # Pipe start (outer boundayr conditions)
        if i == 0:
            V2 = V0[i+1]; H2 = H0[i+1]
            if utype[0] == 'Reservoir' or  utype[0] == 'Tank':
                H[i], V[i] = rev_end (H2, V2, H[i], i, a, g, f, D, dt)
            elif utype[0] == 'Valve':
                H[i], V[i] = valve_end (H2, V2, V[i], i, a, g, f, D, dt)
            elif utype[0] == 'Junction':
                H[i], V[i] = dead_end (linkp , H2, V2, i, a, g, f, D, dt)
            elif utype[0] == 'Pump':  #source pump 
                pump[0][1] = calc_parabola_vertex(pump[0][1])
                H[i], V[i] = source_pump(pump[0], linkp, H2, V2, dt, g, [-1])

        # Pipe end  (inner boundary conditions)
        if i == n: 
            V1 = V0[i-1]; H1 = H0[i-1]     # upstream node 
            V2 = V20;     H2 = H20         # downstream nodes
            if dtype[0] == 'Pipe':
                if leak_loc != None and int(p[pn].end_node.id) in leak_loc :
                    leak_A += linkp.end_demand_coeff
                    H[i], V[i] = add_leakage(leak_A, linkp, link2,
                     H1, V1, H2, V2, dt, g, i, [1], np.sign(links2))
                elif burst_loc != None and int(p[pn].end_node.id) in burst_loc :
                    burst_A += linkp.end_demand_coeff
                    H[i], V[i] = add_leakage(burst_A, linkp, link2,
                     H1, V1, H2, V2, dt, g, i, [1], np.sign(links2))
                else :
                    H[i], V[i] = add_leakage(linkp.end_demand_coeff, linkp, link2,
                     H1, V1, H2, V2, dt, g, i, [1], np.sign(links2))

            elif dtype[0] == 'Pump':
                pumpc = calc_parabola_vertex(pump[1])
                H[i], V[i] = pump_node(pumpc, linkp, link2,
                     H1, V1, H2, V2, dt, g, i, [1], np.sign(links2))
                
            elif dtype[0] == 'Valve':
                valvec = valve[1]
                H[i], V[i] = valve_node(valvec, linkp, link2,
                     H1, V1, H2, V2, dt, g, i, [1], np.sign(links2))
                
        # Interior points
        if (i > 0) and (i < n):
            V1 = V0[i-1]; H1 = H0[i-1]
            V2 = V0[i+1]; H2 = H0[i+1]
            
            H[i], V[i] = inner_node(linkp, linkp, 0,
             H1, V1, H2, V2, dt, g, i,[1],[-1])
            
    return H, V 

def right_boundary(linkp, pn, H, V, links1, p, pump, valve, n, dt,
                H0, V0, H10, V10, utype, dtype, 
                leak_loc, leak_A, burst_loc, burst_A):
    """MOC solution for an indivial right boundary pipe.
    
    Parameters
    ----------
    linkp : object
        Current pipe object
    pn : int
        Current pipe ID
    H : numpy.ndarray 
        Head of current pipe at current time step [m]
    V : numpy.ndarray
        Velocity of current pipe at current time step [m/s]
    links1 : list
        Upstream ajacent pipes 
    p : list
        pipe list 
    pump : list
        Characteristics of the pump 
    valve : list
        Characteristics of the valve
    n : int
        Number of discretization of current pipe
    dt : float 
        Time step 
    H0 : numpy.ndarray
        Head of current pipe at previous time step [m]
    V0 : numpy.ndarray
        Velocity of current pipe at previous time step [m/s]
    H10 : list
        Head of left adjacent nodes at previous time step [m]
    V10 : list
        Velocity of left adjacent nodes at previous time step [m/s]
    utype : list 
        Upstream adjacent link type, and if not pipe, their name
    dtype : list
        Downstream adjacent link type, and if not pipe, their name
    leak_loc : list, optional 
        The list of leakage loction, defined by the name of
        the junction node, by default None
    leak_A : float, optional
        Required if leak_loc is defined
        The leakage coefficient of the leakge, by default None 
        Q_leak = leak_A  [ m^3/s/(m H20)^(1/2)] * \sqrt(H)
    burst_loc : list, optional 
        The list of burst loction,
        defined by the name of the junction node, by default None
    burst_A : float, optional 
        Required if burst_loc is defined
        The burst coefficient of the burst, by default None 
        Q_burst_final = burst_A [ m^3/s/(m H20)^(1/2)] * \sqrt(H)
    
    Returns
    -------
    H : numpy.ndarray
        Head results of the current pipe at current time step. [m]
    V : numpy.ndarray
        Velocity results of the current pipe at current time step. [m/s]
    """

    # Properties of current pipe 
    link1 = [p[abs(i)-1] for i in links1]
    f = linkp.roughness              # unitless
    D = linkp.diameter               # m
    g = 9.8                          # m/s^2
    a = linkp.wavev                  # m/s
    n = int(n)                       # spatial discritization 
    
    for i in range(n+1):
        # Pipe start (inner boundary conditions)
        if i == 0:
            V1 = V10; H1 = H10            # upstream node 
            V2 = V0[i+1]; H2 = H0[i+1]    # downstream node
            if utype[0] == 'Pipe':
                if leak_loc != None and int(p[pn].start_node.id) in leak_loc :
                    leak_A += linkp.start_demand_coeff              
                    H[i], V[i] = add_leakage(leak_A, link1, linkp,
                     H1, V1, H2, V2, dt, g, i, np.sign(links1), [-1])
                    
                elif leak_loc != None and int(p[pn].start_node.id) in burst_loc :
                    burst_A += linkp.start_demand_coeff              
                    H[i], V[i] = add_leakage(burst_A, link1, linkp,
                     H1, V1, H2, V2, dt, g, i, np.sign(links1), [-1])
                    
                else :
                    H[i], V[i] = add_leakage(linkp.start_demand_coeff, link1, linkp,
                     H1, V1, H2, V2, dt, g, i, np.sign(links1), [-1])

            elif utype[0] == 'Pump':
                pumpc = calc_parabola_vertex(pump[0])
                H[i], V[i] = pump_node(pumpc, link1, linkp,
                     H1, V1, H2, V2, dt, g, i,  np.sign(links1), [-1])
            elif utype[0] == 'Valve':
                valvec = valve[0]
                H[i], V[i] = valve_node(valvec, link1, linkp,
                     H1, V1, H2, V2, dt, g, i,  np.sign(links1), [-1])      

        # Pipe end (outer boundary conditions )
        if i == n:
            V1 = V0[i-1]; H1 = H0[i-1]
            if dtype[0] == 'Reservoir' or  dtype[0] == 'Tank':
                H[i], V[i] = rev_end (H1, V1, H[i], i, a, g, f, D, dt)
            if  dtype[0] == 'Valve':
                H[i], V[i] = valve_end (H1, V1, V[i], i, a, g, f, D, dt)
            if dtype[0] == 'Junction':
                H[i], V[i] = dead_end (linkp ,H1, V1, i, a, g, f, D, dt)
                
        # Interior points
        if (i > 0) and (i < n):
            V1 = V0[i-1]; H1 = H0[i-1]
            V2 = V0[i+1]; H2 = H0[i+1]
            
            H[i], V[i] = inner_node(linkp, linkp,0,
             H1, V1, H2, V2, dt, g, i,[1],[-1]) 
                  
    return H, V  