# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 14:05:35 2019

@author: lx2347
"""
from __future__ import print_function
import sys 
import math
import numpy as np
from datetime import datetime 
import wntr
from functools import update_wrapper
import warnings
#from numba import jit
from WMOC.geometry import *
import matplotlib.pyplot as plt


def decorator(d):
    "Make function d a decorator: d wraps a function fn."
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

def wmoc(inp_file, dt, tf, valve_to_close, pump_to_operate,valve_op, pump_op,  
                 pressure_zone_bc, T, leak_loc, leak_A, 
                 burst_loc, burst_A, burst_t, 
                 block_loc, block_A):
    
    "tempory function for test purpose."
    startttime = datetime.now()
    wn, npipe = geometry(inp_file)
    wn, dt, Ndis = discretization(wn, npipe, dt)
    
    print('time step %s s' % dt)  
    tn = int(tf/dt) # Total time steps 
    print ('Total Time Step in this simulation', tn)
    
    # valve closure curve
    vo = valvesetting(dt, tf, valve_op)

    # pump operation curve
    po = pumpsetting(dt, tf, pump_op)

    
    t0 = 0 # initial condition calculated at t0     
    if leak_loc != None :         
        for node in leak_loc:
            leak_node = wn.get_node(node)  
            leak_node.add_leak(wn, area=leak_A/np.sqrt(2*9.8)/1000, discharge_coeff=1,
                               start_time=t0)
        leak_loc = [wn.nodes[i].id-1 for i in leak_loc]
                        
    if block_loc != None :
        block_loc = [wn.nodes[i].id-1 for i in block_loc]

    if burst_loc != None :
        burst_loc = [wn.nodes[i].id-1 for i in burst_loc]
    # initial condition     
    wn, H, V = initialize(wn, t0, Ndis, npipe, tn)  

    # network topology 
    links1, links2, utype, dtype = topology(wn,npipe,pressure_zone_bc)    
    
    # MOC transient simulation 
    wnnodes = wn.nodes 
    wnpipes = wn.pipes 
    wnlinks = wn.links

    H, V, tt = MOC(links1, links2, utype, dtype,
                        wnnodes, wnpipes, wnlinks,
                        H, V, npipe, Ndis, dt, tf, t0, 
                        vo, valve_to_close, po, pump_to_operate, T,
                       leak_loc, leak_A, burst_loc, burst_A, burst_t,
                       block_loc, block_A)
    
    simtime = datetime.now() - startttime
    print('Running Time:', simtime)    
    return  wn, H, V, tt    
    

def initialize(wn, t0, Ndis, npipe, tn):
    # initialize H and V for all computational nodes
    H = [0] * npipe  #pipes 
    V = [0] * npipe   
    g = 9.8
    
    # initialize the results list and set initial conditions   
    for _, pipe in wn.pipes():
        pn = int(pipe.id)-1
        H[pn] = np.zeros((int(Ndis[pn]+1),tn), dtype=np.float64)  # initialize the matrix for heads                        
        V[pn] = np.zeros((int(Ndis[pn]+1),tn), dtype=np.float64)  # initialize the matrix for flow velocity
            
    
    sim = wntr.sim.EpanetSimulator(wn)
    results = sim.run_sim() 
    print ("Initial Condition")
    for _, pipe in wn.pipes():      
        pn = int(pipe.id) -1 
        V[pn][:,0] = np.sign(results.link['flowrate'].loc[t0*3600, pipe.name]) *\
                    results.link['velocity'].loc[t0, pipe.name] *\
                    np.ones((int(Ndis[pn]+1)))
        
        H[pn][:,0] = [results.node['head'].loc[t0, pipe.start_node_name] +\
                      i* (results.node['head'].loc[t0, pipe.end_node_name]-
                       results.node['head'].loc[t0, pipe.start_node_name])/int(Ndis[pn]+1) 
                    for i in range(int(Ndis[pn]+1))]

       # calculate demand coefficient        
        pipe.start_demand_coeff = lambda: None
        pipe.end_demand_coeff = lambda: None 
        try:
            start_demand_coeff = wn.nodes[pipe.start_node_name].demand_timeseries_list.at(t0)/ np.sqrt(H[pn][0,0])
        except :
            start_demand_coeff = 0. 
            
        try: 
            end_demand_coeff =wn.nodes[pipe.end_node_name].demand_timeseries_list.at(t0)/ np.sqrt(H[pn][-1,0])
        except :
            end_demand_coeff = 0. 
            
        setattr(pipe, 'start_demand_coeff',start_demand_coeff*1000 )
        setattr(pipe, 'end_demand_coeff',end_demand_coeff*1000 )        
        
        # tolerance for velocity and head change
        if abs(V[pn][0,0]) >= 1e-5 and abs(H[pn][0,0] - H[pn][-1,0]) >=1e-4:
            pipe.roughness = abs(H[pn][0,0]-H[pn][-1,0]) / (pipe.length/pipe.diameter)/\
                        (V[pn][0,0]**2/2/g)
        else:
            pipe.roughness = 0 
            
        if pipe.roughness >0.08:
            warnings.warn("%s :the friction coefficient %s is too large. \
                          The D-W coeff has been set to 0.03 "%(pipe.name, pipe.roughness))
            pipe.roughness = 0.03
            
#        print ('Pipe', pipe.name, pipe.start_demand_coeff, pipe.end_demand_coeff)
        
    return wn,H, V




def valve_curve(s,valve='Gate'):
    percent_open = np.linspace(100,0,11)
    # loss coefficients for a gate valve
    kl = [1/0.2, 2.50, 1.25, 0.625, 0.333, 0.17,
          0.100, 0.0556, 0.0313, 0.0167, 0.0]
    k = np.interp(s,  percent_open[::-1], kl[::-1])
    return k

@memo 
def MOC(links1, links2, utype, dtype, wnnodes, wnpipes, wnlinks,
        H, V, npipe, Ndis, dt, tf, t0, vo, valve_to_close, po, pump_to_operate,
        T, leak_loc, leak_A, burst_loc, final_burst, burst_t, block_loc, block_A):
    """ Apply Method of Charateristics to the given network, 
    and return results of pressure head and flowrate at each 
    computation nodes.
    Output:
        H[npipe][nnode,tn]"""
          
        
    tt = ['x']
    tt.append(0)

    tn = int(tf/dt) # Total time steps 
    # determine which node of the adjacant pipe should be call:
    # if the adjacant pipe is entering the junction, then -2
    # if the adjacent pipe is leaving the junction, then 1
    a = {1:-2, -1:1}    
    # generat a list of pipe 
    p = []
    for _, pipe in wnpipes():
        p.append(pipe)
    """ Start Claculation"""      
    for ts in range(1,tn):       
        t = ts*dt
#        print ("Calculation time", t)
        tt.append(t) 
        
        # calculate the burst area in the current time
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
        for _,pipe in wnpipes():
            pn = pipe.id-1         
#            print ('Calculating pipe #',pipe.name)
            """Assumption: when a pipe is connected with a pump or valve, 
                the connection is not branch junction.
            """
            # inner pipes 
            if links1[pn] and links2[pn] and \
                links1[pn] != ['End'] and links2[pn] != ['End']:
                
                pump = [0,0]; valve = [0,0]
                if utype[pn][0] == 'Pump': 
                    pump[0] = wnlinks[utype[pn][1]].get_pump_curve().points
                    if utype[pn][1] in pump_to_operate:
                        pump[0]=[(i*po[ts],j*po[ts]**2) for (i,j) in pump[0]]
                        
                elif utype[pn][0] == 'Valve':                                         
                    if utype[pn][1] in valve_to_close:
                        valve[0] = valve_curve(vo[ts]*100)
                    else :  
                         valve[0] = valve_curve(100)
                if dtype[pn][0] == 'Pump':                    
                    pump[1] = wnlinks[dtype[pn][1]].get_pump_curve().points                
                    if dtype[pn][1] in pump_to_operate:
                        pump[1]=[(i*po[ts],j*po[ts]**2) for (i,j) in pump[1]]
   
                elif dtype[pn][0] == 'Valve':
                    if dtype[pn][1] in valve_to_close:
                        valve[1] = valve_curve(vo[ts]*100)
                    else :                        
                         valve[1] = valve_curve(100)
                # demand at the start node and end node
#                demand = [wnnodes[pipe.start_node_name].demand_timeseries_list.at(t0+t),
#                      wnnodes[pipe.end_node_name].demand_timeseries_list.at(t0+t)]

                HP, VP = inner_pipe(pipe, pn,
                     H[pn][:,ts], V[pn][:,ts],
                     links1[pn], links2[pn], utype[pn], dtype[pn], p,
                     np.asscalar(Ndis[pn]), dt,
                     H[pn][:,ts-1], V[pn][:,ts-1], 
                     [H[abs(i)-1][a[np.sign(i)],ts-1] for i in links1[pn]],
                     [V[abs(i)-1][a[np.sign(i)],ts-1] for i in links1[pn]], 
                     [H[abs(i)-1][a[np.sign(i)],ts-1] for i in links2[pn]],
                     [V[abs(i)-1][a[np.sign(i)],ts-1] for i in links2[pn]],
                     pump, valve, 
                     leak_loc, leak_A, burst_loc, burst_A, 
                     block_loc, block_A)    
                H[pn][:,ts], V[pn][:,ts] = HP, VP
                                
            # left boundary pipe           
            elif not links1[pn] or links1[pn] == ['End']:
                pump = [0,0]; valve = [0,0]
#                 LEFT BOUNDARY 
                if utype[pn][0] == 'Reservoir':                    
                    H[pn][0,:]   =  wnnodes[utype[pn][1]].base_head # head of reservoir 
                elif utype[pn][0] == 'Tank':
                    H[pn][0,:]   =  wnnodes[utype[pn][1]].head # head of tank 
                elif utype[pn][0] == 'Junction':
                    V[pn][0,ts] = V[pn][0,0]   
                elif utype[pn][0] == 'PressureZone':
                    V[pn][0,ts] = V[pn][0,0] 
                    T =T 
                elif utype[pn][0] == 'Valve':
                    if utype[pn][1] in valve_to_close:
                        V[pn][0,ts]   = V[pn][0,0] *vo[ts]  # valve velocity condition 
                    else :
                        V[pn][0,:]   = V[pn][0,0]
                elif utype[pn][0] == 'Pump': #source pump         
                    pump[0] = [wnlinks[utype[pn][1]].start_node.base_head,
                         wnlinks[utype[pn][1]].get_pump_curve().points]
                      
                    if utype[pn][1] in pump_to_operate:
                        pump[0][1]=[(i*po[ts],j*po[ts]**2) for (i,j) in pump[0][1]]
                        print('Operating pump %s' %utype[pn][1])
                else:
                     warnings.warn ('Pipe %s miss %s upstream.' %(pipe, utype[pn][0]))   
                    
                # RIGHT BOUNDARY    
                if dtype[pn][0] == 'Pump':                    
                    pump[1] = wnlinks[dtype[pn][1]].get_pump_curve().points                
                    if dtype[pn][1] in pump_to_operate:
                        pump[1]=[(i*po[ts],j*po[ts]**2) for (i,j) in pump[1]]
                        
                elif dtype[pn][0] == 'Valve':
                    if dtype[pn][1] in valve_to_close:
                        valve[1] = valve_curve(vo[ts]*100)
                    else :                        
                         valve[1] = valve_curve(100)
                                           
                # demand at the end node
#                demand = wnnodes[pipe.end_node_name].demand_timeseries_list.at(t0+t)
                HP, VP = left_boundary(pipe, pn, 
                     H[pn][:,ts], V[pn][:,ts],
                     links2[pn], p, T, pump, valve,
                     np.asscalar(Ndis[pn]), dt,
                     H[pn][:,ts-1], V[pn][:,ts-1], 
                     [H[abs(i)-1][a[np.sign(i)],ts-1] for i in links2[pn]],
                     [V[abs(i)-1][a[np.sign(i)],ts-1] for i in links2[pn]],
                     utype[pn], dtype[pn], leak_loc, leak_A, burst_loc, burst_A, 
                     block_loc, block_A) 
                H[pn][:,ts], V[pn][:,ts] = HP, VP
               
            #  right boundary pipe
            elif not links2[pn] or links2[pn] == ['End']:
                pump = [0,0]; valve = [0,0] 
                # RIGHT boundary
                if dtype[pn][0] == 'Reservoir':
                    H[pn][-1,:]   =  wnnodes[dtype[pn][1]].base_head # head of reservoir 
                elif dtype[pn][0] == 'Tank':
                    H[pn][-1,:]   =  wnnodes[dtype[pn][1]].head # head of tank 
                elif dtype[pn][0] == 'Junction':
                    V[pn][-1,ts] = V[pn][-1,0] 
                elif dtype[pn][0] == 'PressureZone':
                    V[pn][-1,ts] = V[pn][-1,0] 
                    T = T 
                elif dtype[pn][0] == 'Valve':
                    if dtype[pn][1] in valve_to_close:
                        V[pn][-1,ts] = V[pn][-1,0]*vo[ts]  # valve velocity condition 
                    else :
                        V[pn][-1,:] = V[pn][-1,0]
                else :
                     warnings.warn('Pipe %s miss %s downstream.' %(pipe, dtype[pn][0]))
                # LEFT boundary   
                if utype[pn][0] == 'Pump': 
                    pump[0] = wnlinks[utype[pn][1]].get_pump_curve().points
                    if utype[pn][1] in pump_to_operate:
                        pump[0]=[(i*po[ts],j*po[ts]**2) for (i,j) in pump[0]]

                        
                elif utype[pn][0] == 'Valve':                                         
                    if utype[pn][1] in valve_to_close:
                        valve[0] = valve_curve(vo[ts]*100)
                    else :  
                         valve[0] = valve_curve(100)
                # demand at the end node
#                demand = wnnodes[pipe.start_node_name].demand_timeseries_list.at(t0+t) 

                HP, VP = right_boundary(pipe, pn,  
                     H[pn][:,ts], V[pn][:,ts],
                     links1[pn], p, T, pump, valve, 
                     np.asscalar(Ndis[pn]), dt,
                     H[pn][:,ts-1], V[pn][:,ts-1], 
                     [H[abs(i)-1][a[np.sign(i)],ts-1] for i in links1[pn]],
                     [V[abs(i)-1][a[np.sign(i)],ts-1] for i in links1[pn]],                
                     utype[pn], dtype[pn], leak_loc, leak_A, burst_loc, burst_A, 
                     block_loc, block_A)   
                H[pn][:,ts], V[pn][:,ts] = HP, VP
                
    return H, V, tt[1:]

def inner_pipe (linkp, pn,  H, V, links1, links2, utype, dtype, p, n, dt, 
                H0, V0, H10, V10, H20, V20, pump, valve,
                leak_loc, leak_A, burst_loc, burst_A, block_loc, block_A):
    """  Transient flow analysis in a single pipe using 
    Method of Characteristics (MOC) 
    Input:
        results from last time step at current pipe (H0,V0),
        L1 (H10,V10), R2 (H20,V20)
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
                    H[i], V[i] = add_leakage(leak_A, link1, linkp,
                     H1, V1, H2, V2, dt, g, i,  np.sign(links1), [-1])
                    
                elif burst_loc != None and int(p[pn].start_node.id) in burst_loc:
                    H[i], V[i] = add_leakage(burst_A, link1, linkp,
                     H1, V1, H2, V2, dt, g, i,  np.sign(links1), [-1])
                    
                elif block_loc != None and int(p[pn].start_node.id) in block_loc:
                    H[i], V[i] = add_blockage(leak_A, V1, H1, link1, 
                     V2, H2, linkp, i, g, dt)
                else :
                    H[i], V[i] = add_leakage(linkp.start_demand_coeff, link1, linkp,
                     H1, V1, H2, V2, dt, g, i,  np.sign(links1), [-1])
#                    H[i], V[i] = inner_boundary(link1, linkp, demand[0],
#                     H1, V1, H2, V2, dt, g, i,  np.sign(links1), [-1])
            elif utype[0] == 'Pump':
                pumpc = calc_parabola_vertex(pump[0])
                H[i], V[i] = pump_boundary(pumpc, link1, linkp,
                     H1, V1, H2, V2, dt, g, i,  np.sign(links1), [-1])
            elif utype[0] == 'Valve':
                valvec = valve[0]
                H[i], V[i] = valve_boundary(valvec, link1, linkp,
                     H1, V1, H2, V2, dt, g, i,  np.sign(links1), [-1])
        # Pipe end  
        if i == n:
            V1 = V0[i-1];    H1 = H0[i-1]
            V2 = V20;        H2 = H20  
            if dtype[0] == 'Pipe':
                if leak_loc != None and int(p[pn].end_node.id) in leak_loc:
                    H[i], V[i] = add_leakage(leak_A,linkp, link2,
                     H1, V1, H2, V2, dt, g, i, [1], np.sign(links2))
                elif burst_loc != None and int(p[pn].end_node.id) in burst_loc:
                    H[i], V[i] = add_leakage(burst_A,linkp, link2,
                     H1, V1, H2, V2, dt, g, i, [1], np.sign(links2))
                elif block_loc != None and int(p[pn].end_node.id) in block_loc:
                    H[i], V[i] = add_blockage(block_A, V1, H1, linkp, 
                     V2, H2, link2, i, g, dt)
                else :
                    H[i], V[i] = add_leakage(linkp.end_demand_coeff,linkp, link2,
                     H1, V1, H2, V2, dt, g, i, [1], np.sign(links2))
                    
#                    H[i], V[i] = inner_boundary(linkp, link2, demand[1],
#                     H1, V1, H2, V2, dt, g, i, [1], np.sign(links2)) 
                    
            elif dtype[0] == 'Pump':
                pumpc = calc_parabola_vertex(pump[1])
                H[i], V[i] = pump_boundary(pumpc, linkp, link2,
                     H1, V1, H2, V2, dt, g, i, [1], np.sign(links2))
                
            elif dtype[0] == 'Valve':
                valvec = valve[1]
                H[i], V[i] = valve_boundary(valvec, linkp, link2,
                     H1, V1, H2, V2, dt, g, i, [1], np.sign(links2))
                
        # Interior points
        if (i > 0) and (i < n):
            V1 = V0[i-1]; H1 = H0[i-1]
            V2 = V0[i+1]; H2 = H0[i+1]
                
            H[i], V[i] = inner_boundary(linkp, linkp, 0,
             H1, V1, H2, V2, dt, g, i,[1],[-1])
                  
    return H, V 


def left_boundary(linkp, pn, H, V, links2, p, T, pump, valve, n, dt, 
                  H0, V0, H20, V20, utype, dtype,
                  leak_loc, leak_A, burst_loc, burst_A, 
                  block_loc, block_A ) :
    """ Left boundary pipe"""
    
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
            elif  utype[0] == 'PressureZone':
                H[i], V[i] = pressure_zone (H2, V2, H0[i], T, i, a, g, f, D, dt)
        # Pipe end  (inner boundary conditions)
        if i == n: 
            V1 = V0[i-1]; H1 = H0[i-1]     # upstream node 
            V2 = V20;     H2 = H20         # downstream nodes
            if dtype[0] == 'Pipe':
                if leak_loc != None and int(p[pn].end_node.id) in leak_loc :
                    H[i], V[i] = add_leakage(leak_A, linkp, link2,
                     H1, V1, H2, V2, dt, g, i, [1], np.sign(links2))
                elif burst_loc != None and int(p[pn].end_node.id) in burst_loc :
                    H[i], V[i] = add_leakage(burst_A, linkp, link2,
                     H1, V1, H2, V2, dt, g, i, [1], np.sign(links2))
                elif block_loc != None and int(p[pn].end_node.id) in block_loc :
                    H[i], V[i] = add_blockage(block_A, V1, H1, linkp, 
                     V2, H2, link2, i, g, dt)
                else :
                    H[i], V[i] = add_leakage(linkp.end_demand_coeff, linkp, link2,
                     H1, V1, H2, V2, dt, g, i, [1], np.sign(links2))
#                    H[i], V[i] = inner_boundary(linkp, link2, demand, 
#                     H1, V1, H2, V2, dt, g, i, [1], np.sign(links2))
            elif dtype[0] == 'Pump':
                pumpc = calc_parabola_vertex(pump[1])
                H[i], V[i] = pump_boundary(pumpc, linkp, link2,
                     H1, V1, H2, V2, dt, g, i, [1], np.sign(links2))
                
            elif dtype[0] == 'Valve':
                valvec = valve[1]
                H[i], V[i] = valve_boundary(valvec, linkp, link2,
                     H1, V1, H2, V2, dt, g, i, [1], np.sign(links2))
                
        # Interior points
        if (i > 0) and (i < n):
            V1 = V0[i-1]; H1 = H0[i-1]
            V2 = V0[i+1]; H2 = H0[i+1]
            
            H[i], V[i] = inner_boundary(linkp, linkp, 0,
             H1, V1, H2, V2, dt, g, i,[1],[-1])
            
    return H, V 

def right_boundary(linkp, pn, H, V, links1, p, T, pump, valve, n, dt,
                H0, V0, H10, V10, utype, dtype, 
                leak_loc, leak_A, burst_loc, burst_A, 
                block_loc, block_A):
    """ Right boundary pipe """
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
                    H[i], V[i] = add_leakage(leak_A, link1, linkp,
                     H1, V1, H2, V2, dt, g, i, np.sign(links1), [-1])
                    
                elif leak_loc != None and int(p[pn].start_node.id) in burst_loc :              
                    H[i], V[i] = add_leakage(burst_A, link1, linkp,
                     H1, V1, H2, V2, dt, g, i, np.sign(links1), [-1])
                    
                elif block_loc != None and int(p[pn].start_node.id) in block_loc :                
                    H[i], V[i] = add_blockage(block_A, V1, H1, link1, 
                     V2, H2, linkp, i, g, dt)
                else :
                    H[i], V[i] = add_leakage(linkp.start_demand_coeff, link1, linkp,
                     H1, V1, H2, V2, dt, g, i, np.sign(links1), [-1])
#                    H[i], V[i] = inner_boundary(link1, linkp, demand,
#                     H1, V1, H2, V2, dt, g, i, np.sign(links1), [-1])
            elif utype[0] == 'Pump':
                pumpc = calc_parabola_vertex(pump[0])
                H[i], V[i] = pump_boundary(pumpc, link1, linkp,
                     H1, V1, H2, V2, dt, g, i,  np.sign(links1), [-1])
            elif utype[0] == 'Valve':
                valvec = valve[0]
                H[i], V[i] = valve_boundary(valvec, link1, linkp,
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
            if dtype[0] == 'PressureZone':
                H[i], V[i] = pressure_zone (H1, V1, H0[i], T, i, a, g, f, D, dt)
                
        # Interior points
        if (i > 0) and (i < n):
            V1 = V0[i-1]; H1 = H0[i-1]
            V2 = V0[i+1]; H2 = H0[i+1]
            
            H[i], V[i] = inner_boundary(linkp, linkp,0,
             H1, V1, H2, V2, dt, g, i,[1],[-1]) 
                  
    return H, V  

def calc_parabola_vertex(points):
    '''Adapted and modifed to get the unknowns for defining a parabola:'''
    [(x1,y1),(x2,y2),(x3,y3)] = points
    denom = (x1-x2) * (x1-x3) * (x2-x3)
    A     = (x3 * (y2-y1) + x2 * (y1-y3) + x1 * (y3-y2)) / denom
    B     = (x3*x3 * (y1-y2) + x2*x2 * (y3-y1) + x1*x1 * (y2-y3)) / denom
    C     = (x2 * x3 * (x2-x3) * y1+x3 * x1 * (x3-x1) * y2+x1 * x2 * (x1-x2) * y3) / denom
    return A,B,C



    
def inner_boundary(link1, link2, demand, H1, V1, H2, V2, dt, g, nn, s1, s2):
    """Deal with inner nodes.
    Coded according to Larock, Jeppson, and Watters's
    Hydraulics of Pipeline Systems"""  
    try :
        list(link1) 
    except:
        link1 = [link1]
        V1 = [V1] ; H1 = [H1]
            
    try :
        list(link2)
    except:
        link2 = [link2]  
        V2 = [V2] ; H2 = [H2]

    # property of left adjacent pipe   
    f1 = [link1[i].roughness  for i in range(len(link1))]   # unitless
    D1 = [link1[i].diameter  for i in range(len(link1))]    # m 
    a1 = [link1[i].wavev  for i in range(len(link1))] # m/s 
    A1 = [math.pi * D1[i]**2. / 4.  for i in range(len(link1))]   #m^2
    C1 = np.zeros((len(link1),2), dtype=np.float64)   
    theta1 = [link1[i].theta for i in range((len(link1)))]
    

    for i in range(len(link1)):
        C1[i,0] = s1[i]*V1[i] + g/a1[i]*H1[i] - (s1[i]*f1[i]*dt
          /2./D1[i]*V1[i]*abs(V1[i])) + g/a1[i]* dt *V1[i]*theta1[i]
        C1[i,1] = g/a1[i]
        
#    H-W coefficients given  
#    C1[i,0] = s1[i]*V1[i] + g/a1[i]*H1[i] - (s1[i]*6.8405*dt*g
#      /f1[i]**1.852/D1[i]**1.166*V1[i]*abs(V1[i])**0.852) +g/a1[i]* dt *V1[i]*theta1[i]       
#    C1[i,1] = g/a1[i]
        
    # property of right adjacent pipe    
    f2 = [link2[i].roughness  for i in range(len(link2))]      # unitless
    D2 = [link2[i].diameter  for i in range(len(link2))]       #m  
    a2 = [link2[i].wavev  for i in range(len(link2))] #m/s
    A2 = [math.pi * D2[i]**2. / 4.  for i in range(len(link2))]    #m^2
    C2 = np.zeros((len(link2),2),dtype=np.float64)
    theta2 = [link2[i].theta for i in range((len(link2)))]
        
    for i in range(len(link2)):
        C2[i,0] = s2[i]*V2[i] + g/a2[i]*H2[i] - (s2[i]*f2[i]*
          dt/2./D2[i]*V2[i]*abs(V2[i])) + g/a2[i]* dt *V2[i]*theta2[i]
        C2[i,1] = g/a2[i]
    
#    H-W coefficients given
#    C2[i,0] = s2[i]*V2[i] + g/a2[i]*H2[i] - (s2[i]*6.8405*dt*g
#      /f2[i]**1.852/D2[i]**1.166*V2[i]*abs(V2[i])**0.852) + g/a2[i]* dt *V2[i]*theta2[i]
#    C2[i,1] = g/a2[i]
            
    if link1 == link2 : # inner node of one pipe 
        HP = ((np.dot(C1[:,0], A1) + np.dot(C2[:,0],A2)) / 
         (np.dot(C1[:,1], A1) + np.dot(C2[:,1],A2)))

    else : # junction 
        HP = (((np.dot(C1[:,0], A1) + np.dot(C2[:,0],A2)) - demand)/ 
         (np.dot(C1[:,1], A1) + np.dot(C2[:,1],A2)))        
    if nn == 0:  # pipe start 
        VP = np.float64(-C2[:,0]+ C2[:,1]*HP)
    else:        # pipe end 
        VP = np.float64(C1[:,0] - C1[:,1]*HP)
    return HP, VP 


def valve_boundary(KL_inv, link1, link2, H1, V1, H2, V2, dt, g, nn, s1, s2):
#    print (KL_inv)
    """Deal with inline valve.
    Coded according to Larock, Jeppson, and Watters's
    Hydraulics of Pipeline Systems"""  
    try :
        list(link1) 
    except:
        link1 = [link1]
        V1 = [V1] ; H1 = [H1]
            
    try :
        list(link2)
    except:
        link2 = [link2]  
        V2 = [V2] ; H2 = [H2]   
    # property of left adjacent pipe   
    f1 = [link1[i].roughness  for i in range(len(link1))]   # unitless
    D1 = [link1[i].diameter  for i in range(len(link1))]    # m 
    a1 = [link1[i].wavev  for i in range(len(link1))] # m/s 
    A1 = [math.pi * D1[i]**2. / 4.  for i in range(len(link1))]   #m^2
    C1 = np.zeros((len(link1),2), dtype=np.float64)   
    theta1 = [link1[i].theta for i in range((len(link1)))]

    for i in range(len(link1)):
        C1[i,0] = s1[i]*V1[i] + g/a1[i]*H1[i] - (s1[i]*f1[i]*dt
          /2./D1[i]*V1[i]*abs(V1[i])) + g/a1[i]* dt *V2[i]*theta1[i]
        C1[i,1] = g/a1[i]
# H-W coefficients given  
#    for i in range(len(link1)):
#        C1[i,0] = s1[i]*V1[i] + g/a1[i]*H1[i] - (s1[i]*6.8405*dt*g
#          /f1[i]**1.852/D1[i]**1.166*V1[i]*abs(V1[i])**0.852) + g/a1[i]* dt *V2[i]*theta1[i]       
#        C1[i,1] = g/a1[i]

        
    # property of right adjacent pipe    
    f2 = [link2[i].roughness  for i in range(len(link2))]      # unitless
    D2 = [link2[i].diameter  for i in range(len(link2))]       #m  
    a2 = [link2[i].wavev  for i in range(len(link2))]          #m/s
    A2 = [math.pi * D2[i]**2. / 4.  for i in range(len(link2))]    #m^2
    C2 = np.zeros((len(link2),2),dtype=np.float64)
    theta2 = [link2[i].theta for i in range((len(link2)))]

    for i in range(len(link2)):    
        C2[i,0] = s2[i]*V2[i] + g/a2[i]*H2[i] - (s2[i]*f2[i]*
          dt/2./D2[i]*V2[i]*abs(V2[i])) + g/a2[i]*dt*V2[i]*theta2[i]
        C2[i,1] = g/a2[i]
#   H-W coefficients given
#    for i in range(len(link2)):
#        C2[i,0] = s2[i]*V2[i] + g/a2[i]*H2[i] - (s2[i]*6.8405*dt*g
#          /f2[i]**1.852/D2[i]**1.166*V2[i]*abs(V2[i])**0.852)+ g/a2[i]*dt*V2[i]*theta2[i]
#        C2[i,1] = g/a2[i]
   
    # parameters of the quatratic polinomial
    aq = 1
    bq = 2*g*KL_inv* (A2[0]/A1[0]/C1[0,1] + 1/C2[0,1])
    cq = 2*g*KL_inv* (C2[0,0]/C2[0,1] - C1[0,0]/C1[0,1])

    # solve the quadratic equation
    delta = bq**2 - 4*aq*cq  

    if delta >= 0:
        VP = (-bq + np.sqrt(delta))/(2*aq)
    elif delta > -1.0e-7 and delta <0 :
        VP = (-bq)/(2*aq)
    else:
        VP = (-bq)/(2*aq)
        sys.exit('Error: The quadratic equation has no real solution (valve)')
    
    if VP >=0 : # positive flow            
        if nn == 0:  # pipe start       
            VP = VP
            HP = (C2[0,0] + VP) / C2[0,1]
        else:        # pipe end 
            VP = VP*A2[0]/A1[0]
            HP = (C1[0,0] - VP) / C1[0,1]      
            
    else : # reverse flow
        # reconsruct the quadratic equation
        # parameters of the quatratic polinomial
        aq = 1
        bq = 2*g*KL_inv* (-A1[0]/A2[0]/C2[0,1]-1/C1[0,1])
        cq = 2*g*KL_inv* (-C2[0,0]/C2[0,1]+C1[0,0]/C1[0,1])
        
        # solve the quadratic equation
        delta = bq**2 - 4*aq*cq  
        
        if delta >= 0:
            VP = (-bq - np.sqrt(delta))/(2*aq)
        elif delta > -1.0e-7 and delta <0 :
            VP = (-bq)/(2*aq)
        else:
            VP = (-bq)/(2*aq)
            sys.exit('Error: The quadratic equation has no real solution (valve)')
        
        if nn == 0:  # pipe start       
            VP = VP*A1[0]/A2[0]
            HP = (C2[0,0] + VP ) / C2[0,1]
        else:        # pipe end 
            VP = VP
            HP = (C1[0,0] - VP) / C1[0,1]      
    return HP, VP 
 
    
def pump_boundary(pumpc,link1, link2, H1, V1, H2, V2, dt, g, nn, s1, s2):
#    print ('s', s1, s2 )
    """Deal with pump boundary conditions.
    Coded according to Larock, Jeppson, and Watters's
    Hydraulics of Pipeline Systems"""  

    try :
        list(link1) 
    except:
        link1 = [link1]
        V1 = [V1] ; H1 = [H1]
            
    try :
        list(link2)
    except:
        link2 = [link2]  
        V2 = [V2] ; H2 = [H2]
    
    # property of left adjacent pipe   
    f1 = [link1[i].roughness  for i in range(len(link1))]       # unitless
    D1 = [link1[i].diameter  for i in range(len(link1))]        # m 
    a1 = [link1[i].wavev  for i in range(len(link1))]           # m/s 
    A1 = [math.pi * D1[i]**2. / 4.  for i in range(len(link1))] # m^2
    C1 = np.zeros((len(link1),2), dtype=np.float64)   
    theta1 = [link1[i].theta for i in range((len(link1)))]
 
    # D-W coefficients given 
    for i in range(len(link1)):
        C1[i,0] = s1[i]*V1[i] + g/a1[i]*H1[i] - (s1[i]*f1[i]*dt
          /2./D1[i]*V1[i]*abs(V1[i])) + g/a1[i]*dt*V1[i]*theta1[i]
        C1[i,1] = g/a1[i]
    # H-W coefficients given  
#     for i in range(len(link1)):
#        C1[i,0] = s1[i]*V1[i] + g/a1[i]*H1[i] - (s1[i]*6.8405*dt*g
#          /f1[i]**1.852/D1[i]**1.166*V1[i]*abs(V1[i])**0.852) + g/a1[i]*dt*V1[i]*theta1[i]       
#        C1[i,1] = g/a1[i]
        
    # property of right adjacent pipe    
    f2 = [link2[i].roughness  for i in range(len(link2))]      # unitless
    D2 = [link2[i].diameter  for i in range(len(link2))]       # m  
    a2 = [link2[i].wavev  for i in range(len(link2))]          # m/s
    A2 = [math.pi * D2[i]**2. / 4.  for i in range(len(link2))]# m^2
    C2 = np.zeros((len(link2),2),dtype=np.float64)
    theta2 = [link2[i].theta for i in range(len(link2))]

    for i in range(len(link2)):    
        C2[i,0] = s2[i]*V2[i] + g/a2[i]*H2[i] - (s2[i]*f2[i]*
          dt/2./D2[i]*V2[i]*abs(V2[i])) + g/a2[i]*dt*V2[i]*theta2[i]
        C2[i,1] = g/a2[i]
    # H-W coefficients given
#    for i in range(len(link2)):
#        C2[i,0] = s2[i]*V2[i] + g/a2[i]*H2[i] - (s2[i]*6.8405*dt*g
#          /f2[i]**1.852/D2[i]**1.166*V2[i]*abs(V2[i])**0.852)+ g/a2[i]*dt*V2[i]*theta2[i]
#        C2[i,1] = g/a2[i]
            
    # pump power function 
    ap, bp, cp = pumpc
    ap = ap * A1[0]**2.
    bp = bp * A1[0]
    
    # parameters of the quatratic polinomial
    aq = 1
    bq = 1/ap * (bp - 1/C1[0,1] - A1[0]/C2[0,1]/A2[0])
    cq = 1/ap * (-C2[0,0]/C2[0,1] + C1[0,0]/C1[0,1] + cp)
    
    # solve the quadratic equation
    delta = bq**2. - 4.*aq*cq 
    if delta >= 0:
        VP = (-bq + np.sqrt(delta))/(2*aq)        
    elif delta > -1.0e-7 and delta <0 :
        VP = (-bq)/(2*aq)
    else:
        VP = (-bq)/(2*aq)
        sys.exit('Error: The quadratic equation has no real solution (pump)')
        
    if VP > 0 : # positive flow 
        if nn == 0:  # pipe start 
            VP = VP*A1[0]/A2[0]
            HP = (C2[0,0] + VP ) / C2[0,1]
        else:        # pipe end 
            VP = VP
            HP = (C1[0,0] - VP) / C1[0,1]
    else :
        warnings.warn( "Reverse flow stopped by check valve!")
        VP = 0 
        if nn == 0:  # pipe start 
            HP = (C2[0,0] + VP ) / C2[0,1]
        else :
            HP = (C1[0,0] - VP) / C1[0,1]             
    return HP, VP 

def source_pump(pump, link2, H2, V2, dt, g, s2):
    """Deal with source pump boundary conditions.
    Coded according to Larock, Jeppson, and Watters's
    Hydraulics of Pipeline Systems"""  
    pumpc = pump[1]
    Hsump = pump[0]
    try :
        list(link2)
    except:
        link2 = [link2]  
        V2 = [V2] ; H2 = [H2]

    # property of right adjacent pipe    
    f2 = [link2[i].roughness  for i in range(len(link2))]      # unitless
    D2 = [link2[i].diameter  for i in range(len(link2))]       #m  
    a2 = [link2[i].wavev  for i in range(len(link2))] #m/s
    A2 = [math.pi * D2[i]**2. / 4.  for i in range(len(link2))]    #m^2
    C2 = np.zeros((len(link2),2),dtype=np.float64)
    theta2 = [link2[i].theta for i in range(len(link2))]

    # D-W coefficients given 
    for i in range(len(link2)):    
        C2[i,0] = s2[i]*V2[i] + g/a2[i]*H2[i] - (s2[i]*f2[i]*
          dt/2./D2[i]*V2[i]*abs(V2[i])) + g/a2[i]*dt*V2[i]*theta2[i]
        C2[i,1] = g/a2[i]
    # H-W coefficients given
#    for i in range(len(link2)):
#        C2[i,0] = s2[i]*V2[i] + g/a2[i]*H2[i] - (s2[i]*6.8405*dt*g
#          /f2[i]**1.852/D2[i]**1.166*V2[i]*abs(V2[i])**0.852)+ g/a2[i]*dt*V2[i]*theta2[i]
#        C2[i,1] = g/a2[i]
            
    # pump power function 
    ap, bp, cp = pumpc
    ap = ap * A2[0]**2.
    bp = bp * A2[0]
    
    # parameters of the quatratic polinomial
    aq = ap * C2[0,1]**2.
    bq = bp*C2[0,1] - 2.*ap*C2[0,0]*C2[0,1] - 1
    cq = ap*C2[0,0]**2. - bp*C2[0,0] + Hsump + cp
    
    # solve the quadratic equation
    delta = bq**2. - 4.*aq*cq 
    if delta >= 0:
        HP = (-bq - np.sqrt(delta))/(2*aq)        
    elif delta > -1.0e-7 and delta <0 :
        HP = (-bq)/(2*aq)
    else:
        HP = (-bq)/(2*aq)
        sys.exit('Error: The quadratic equation has no real solution (pump)')
        
    if HP > Hsump: 
        VP = np.float64(-C2[0,0]+ C2[0,1]*HP)
    else :
        HP = Hsump
        VP = np.float64(-C2[0,0]+ C2[0,1]*HP) 
        
    if VP <= 0 : # positive flow 
        warnings.warn( "Reverse flow stopped by check valve!")
        VP = 0       
        HP = (C2[0,0] + VP ) / C2[0,1]
            
    return HP, VP 



def valve_end(H1, V1, V, nn, a, g, f, D, dt):
    if nn == 0 :
        HP = H1 + a/g*(V - V1) + a/g*f*dt/(2.*D)*V1*abs(V1)
        VP = V
    else :
        HP = H1 - a/g*(V - V1) - a/g*f*dt/(2.*D)*V1*abs(V1)
        VP = V    
    return HP,VP

def dead_end(linkp, H1, V1, nn, a, g, f, D, dt):
    
    A = math.pi/4. * linkp.diameter**2.
    
    if nn == 0:
        k = linkp.start_demand_coeff/1000.
        aq = 1 
        bq = -a/g*k/A
        cq = a/g *V1 - a/g*f*dt/(2.*D)*V1*abs(V1) - H1
        
        # solve the quadratic equation
        delta = bq**2. - 4.*aq*cq 
        if delta >= 0:
            HP = (-bq - np.sqrt(delta))/(2*aq) 
            HP = HP**2.
        elif delta > -1.0e-7 and delta <0 :
            HP = (-bq)/(2*aq)
            HP = HP**2.
        else:
            HP = (-bq)/(2*aq)
            HP = HP**2.
            sys.exit('Error: The quadratic equation has no real solution (dead end)')            
        VP = V1 - g/a *H1 - f*dt/(2.*D)*V1*abs(V1) + g/a*HP 
    else :
        k = linkp.end_demand_coeff/1000.
        aq = 1 
        bq = a/g*k/A
        cq = -a/g *V1 + a/g*f*dt/(2.*D)*V1*abs(V1) - H1        
        # solve the quadratic equation
        delta = bq**2. - 4.*aq*cq 
        if delta >= 0:
            HP = (-bq + np.sqrt(delta))/(2*aq) 
            HP = HP**2.
        elif delta > -1.0e-7 and delta <0 :
            HP = (-bq)/(2*aq)
            HP = HP**2.
        else:
            HP = (-bq)/(2*aq)
            HP = HP**2.
            sys.exit('Error: The quadratic equation has no real solution (dead end)')
        VP = V1 + g/a *H1 - f*dt/(2.*D)*V1*abs(V1) - g/a*HP 
    return HP,VP

def rev_end( H2, V2, H, nn, a, g, f, D, dt):
    if nn == 0 :
        VP = V2 + g/a*(H - H2) - f*dt/(2.*D)*V2*abs(V2)
        HP = H
    else:
        VP = V2 - g/a*(H - H2) - f*dt/(2.*D)*V2*abs(V2)
        HP = H    
    return HP,VP

def tank_end():
    HP = 0
    VP = 0 
    return HP, VP

def pressure_zone(H2, V2, H, T, nn, a, g, f, D, dt):
    """Deal with pressure-zone B.C."""  
    if nn == 0 :
        HP = (1-T)*H + T*H2
        VP = V2 + g/a*(HP - H2) - f*dt/(2.*D)*V2*abs(V2)
    else:
        HP = (1-T)*H + T*H2
        VP = V2 - g/a*(HP - H2) - f*dt/(2.*D)*V2*abs(V2)
        
#    V =  T/ (2-T) * (V1 - V)+ V
#    if nn == 0 :
#        HP = H1 + a/g*(V - V1) + a/g*f*dt/(2.*D)*V1*abs(V1)
#        VP = V
#    else :
#        HP = H1 - a/g*(V - V1) - a/g*f*dt/(2.*D)*V1*abs(V1)
#        VP = V
    return HP, VP   


def add_leakage(emitter_coef, link1, link2, H1, V1, H2, V2, dt, g, nn, s1, s2):
    """Deal with leaking nodes.
    Coded according to Larock, Jeppson, and Watters's
    Hydraulics of Pipeline Systems"""  
    
    emitter_coef = emitter_coef/1000  # m^3/s//(m H2O)^(1/2)
    
    try :
        list(link1) 
    except:
        link1 = [link1]
        V1 = [V1] ; H1 = [H1]
            
    try :
        list(link2)
    except:
        link2 = [link2]  
        V2 = [V2] ; H2 = [H2]

    # property of left adjacent pipe   
    f1 = [link1[i].roughness  for i in range(len(link1))]   # unitless
    D1 = [link1[i].diameter  for i in range(len(link1))]    # m 
    a1 = [link1[i].wavev  for i in range(len(link1))] # m/s 
    A1 = [math.pi * D1[i]**2. / 4.  for i in range(len(link1))]   #m^2
    C1 = np.zeros((len(link1),2), dtype=np.float64)   
    theta1 = [link1[i].theta for i in range(len(link1))]
    # D-W coefficients given 
    for i in range(len(link1)):
        C1[i,0] = s1[i]*V1[i] + g/a1[i]*H1[i] - (s1[i]*f1[i]*dt
          /2./D1[i]*V1[i]*abs(V1[i]))+ g/a1[i]*dt*V1[i]*theta1[i]
        C1[i,1] = g/a1[i]
# H-W coefficients given 
#     for i in range(len(link1)):
#        C1[i,0] = s1[i]*V1[i] + g/a1[i]*H1[i] - (s1[i]*6.8405*dt*g
#          /f1[i]**1.852/D1[i]**1.166*V1[i]*abs(V1[i])**0.852)+ g/a1[i]*dt*V1[i]*theta1[i]        
#        C1[i,1] = g/a1[i]

        
    # property of right adjacent pipe    
    f2 = [link2[i].roughness  for i in range(len(link2))]      # unitless
    D2 = [link2[i].diameter  for i in range(len(link2))]       #m  
    a2 = [link2[i].wavev  for i in range(len(link2))] #m/s
    A2 = [math.pi * D2[i]**2. / 4.  for i in range(len(link2))]    #m^2
    C2 = np.zeros((len(link2),2),dtype=np.float64)
    theta2 = [link2[i].theta for i in range(len(link2))] 
#   D-W coefficients given 
    for i in range(len(link2)):    
        C2[i,0] = s2[i]*V2[i] + g/a2[i]*H2[i] - (s2[i]*f2[i]*
          dt/2./D2[i]*V2[i]*abs(V2[i])) + g/a2[i]*dt*V2[i]*theta2[i]
        C2[i,1] = g/a2[i]
        
#   H-W coefficients given
#    for i in range(len(link2)):
#        C2[i,0] = s2[i]*V2[i] + g/a2[i]*H2[i] - (s2[i]*6.8405*dt*g
#          /f2[i]**1.852/D2[i]**1.166*V2[i]*abs(V2[i])**0.852)+ g/a2[i]*dt*V2[i]*theta2[i]
#        C2[i,1] = g/a2[i]
    
    a = np.dot(C1[:,0], A1) + np.dot(C2[:,0],A2)
    b = np.dot(C1[:,1], A1) + np.dot(C2[:,1],A2)
    # parameters of the quatratic polinomial
    a1 = b**2.
    b1 = -(2.*a*b +emitter_coef**2)
    c1 = a**2.
    
    # solve the quadratic equation
    delta = b1**2 - 4*a1*c1    
    if delta >= 0:
        HP = (-b1 - np.sqrt(delta))/(2*a1)
    elif delta > -1.0e-7 and delta <0 :
        HP = (-b1)/(2*a1)
    else:
        HP = (-b1)/(2*a1)
        sys.exit('Error: The quadratic equation has no real solution (leakage)')
  
    if nn == 0:  # pipe start 
        VP = np.float64(-C2[:,0]+ C2[:,1]*HP)
    else:        # pipe end 
        VP = np.float64(C1[:,0] - C1[:,1]*HP)
    return HP, VP 


def add_blockage(r, V1, H1,link1, V2, H2, link2, i, g, dt):
    """Add blockage (blocakge rate = r) to pipe nodes (not computation nodes). """
    # property of left ajacent pipe
    f1 = link1.roughness              # unitless
    D1 = link1.diameter               # m
    a1 = np.asscalar(link1.wavev)     # m/s
    A1 = math.pi * D1**2. / 4.
    
    # property of right ajacent pipe 
    f2 = link2.roughness              # unitless
    D2 = link2.diameter               # m
    a2 = np.asscalar(link2.wavev)     # m/s
    A2 = math.pi * D2**2. / 4.
    
    C3 = V2 - g/a2*H2 - f2*dt/2./D2*V2*abs(V2) 
    C4 = g/a2
    C1 = V1 + g/a1*H1 - f1*dt/2./D1*V1*abs(V1)
    C2 = g/a1
    
    A1 = math.pi * D1**2. / 4.
    A2 = math.pi * D2**2. / 4.
    
    H = (C1*A1 - C3*A2) / ((1-r)*C4*A2 + C2*A1)

    if i == 0:
        HP = (1-r) * H
        VP = C3 + C4 *HP
    else :
        HP = H
        VP = C1 - C2* HP   
    
    return HP, VP

     


