# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 14:12:09 2019

@author: lx2347
"""
from __future__ import print_function
import math
import numpy as np
import wntr
import matplotlib.pyplot as plt
import networkx as nx

def geometry(inp_file, wavespeed=1200, plot= True):
    """Read geometry and assign parameters
    
    Parameters
    ----------
    inp_file : wntr.network.model.WaterNetworkModel
        .inp file used for EPAnet simulation 
    wavespeed : int, optional
        Wave speed, by default 1200 [m/s]
    plot : bool, optional
        Plot network, by default True

    Returns
    -------
    wn : wntr.network.model.WaterNetworkModel
        Network with updated parameters 
    npipe : integer
        Number of pipes
    """
       
    wn = wntr.network.WaterNetworkModel(inp_file)
    G = wn.get_graph()
    links = wn.links
    npipe = len(links.pipe_names)
    
#    npipe = len(links)
    print ('Number of pipe: %s' %npipe)
    if isinstance(wavespeed, int):
        wavev = wavespeed * np.ones((npipe, 1))
        
    # assign wave speed to each pipes 
    i= 0 
    for _, pipe in wn.pipes():  
        pipe.wavev = lambda: None        
        setattr(pipe, 'wavev', wavev[i])
        i+=1
        
    # calculate the slope for each pipe 
    for _, pipe in wn.pipes():  
        pipe.theta = lambda: None
        try: 
            theta = np.sin(np.arctan(pipe.end_node.elevation -
                pipe.start_node.elevation)/pipe.length)
        except:
            theta = 0.0
        setattr(pipe, 'theta',theta )
        
    # assign ID to each links, start from 1.      
    i =1
    for _, link in wn.links():  
        link.id = lambda: None 
        setattr(link, 'id', i)
        i+=1 
   
    # assign ID to each links, start from 1.      
    i =1
    for _, node in wn.nodes():  
        node.id = lambda: None 
        setattr(node, 'id', i)
        i+=1     ## Graph the network
        
    if plot: 
        # Compute betweenness centrality
        plt.figure(figsize=(4,4), dpi=80, facecolor='w', edgecolor='k') 
        bet_cen = nx.betweenness_centrality(G)
        wntr.graphics.plot_network(wn, node_attribute=bet_cen, 
                                node_size=30, 
                                title='Betweenness Centrality')
    return  wn, npipe

def topology(wn, npipe):
    """Figure out the topology of the network
    
    Parameters
    ----------
    wn : wntr.network.model.WaterNetworkModel
        .inp file used for EPAnet simulation
    npipe : integer
        Number of pipes

    Returns
    -------
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
    """
    G = wn.get_graph()
    
    length = wn.query_link_attribute('length')
    G.weight_graph(link_attribute = length)

    # add 'id' attribute to networkx links 
    i =1
    for ln, link in wn.links():          
        G.edges[link.start_node_name,link.end_node_name, ln]['id'] = i
        i+=1 
        
    # allocate the parameters
    links1 = [0] * len(wn.links) 
    links2 = [0] * len(wn.links)  
    utype = [('Pipe',0)] * npipe
    dtype = [('Pipe',0)] * npipe
 
    # Adjcant pipes for each pipe IN:+; OUT:-   
    for _, link in wn.links():
        pn = link.id
        links1[int(pn)-1] = [int(p['id'])
                           for _, attr in G.pred[link.start_node_name].items() 
                           for _,p in attr.items()
                           if p['id'] != pn] 
        
        for _, attr in G.succ[link.start_node_name].items() :
            for _,p in attr.items():
                if p['id'] != pn:
                    links1[int(pn)-1].append(-1* int(p['id']))
               
                
        # right (end) adjcant pipes 
        links2[int(pn)-1] = [int(p['id']) 
                       for _, attr in G.pred[link.end_node_name].items() 
                       for _,p in attr.items()
                       if p['id'] != pn]  
        
        for _, attr in G.succ[link.end_node_name].items() :
            for _,p in attr.items():
                if p['id'] != pn:
                    links2[int(pn)-1].append(-1*int(p['id']))
                    
    #figure out downstream type and upstream type                
    for _,pipe in wn.pipes():
        pn = pipe.id-1 
        if links1[pn] :   
            if max(map(abs, links1[pn])) > npipe:   
                utype[pn] = [(l.link_type,l.name)
                         for _,l in wn.links() 
                         if l.id == abs(links1[pn][0])][0] 
                
                if links1[abs(links1[pn][0])-1] and links2[abs(links1[pn][0])-1]:
                    links1[pn] = [i 
                    for i in [links1[abs(links1[pn][0])-1], links2[abs(links1[pn][0])-1]]
                    if abs(i[0]) -1 != pn][0]         
                else:
                    links1[pn] = ['End']
                    
        else:
            utype[pn] = (wn.nodes[pipe.start_node_name].node_type, 
                         wn.nodes[pipe.start_node_name])

        if links2[pn] :            
            if max(map(abs, links2[pn])) > npipe:
                dtype[pn] = [(l.link_type,l.name) 
                         for _,l in wn.links() 
                         if l.id == abs(links2[pn][0])][0]   
                 
                if links1[abs(links2[pn][0])-1] and links2[abs(links2[pn][0])-1]:
                    links2[pn] = [i 
                    for i in [links1[abs(links2[pn][0])-1], links2[abs(links2[pn][0])-1]]
                    if abs(i[0]) -1 != pn][0]  
                else:
                    links2[pn] = ['End']
                           
        else:
            dtype[pn] = (wn.nodes[pipe.end_node_name].node_type,
                         wn.nodes[pipe.end_node_name])
        
    return links1, links2, utype, dtype

def max_time_step(wn):
    """Detrtmine the maximum time step based onCourant's criteria.
    
    Parameters
    ----------
    wn : wntr.network.model.WaterNetworkModel
        Network 
    
    Returns
    -------
    max_dt : float 
        Maximum time step allowed for this network
    """
    max_dt = math.inf
 
    for _, pipe in wn.pipes():
        dt = pipe.length / (2. * pipe.wavev)
        if max_dt > dt :
            max_dt = dt           
    return max_dt

def cal_N(wn, npipe, dt):
    """Determine the number of computation uites ($N_i$) for each pipes.
    
    Parameters
    ----------
    wn : wntr.network.model.WaterNetworkModel
        Network
    npipe : integer
        Number of pipes
    dt : float 
        Time step for transient simulation
    """
    N = np.zeros((npipe,1))

    for _, pipe in wn.pipes() :

        N[int(pipe.id)-1] =  int(2*np.int(pipe.length/ (2. * pipe.wavev *dt)))
    return N


def adjust_wavev( wn, N):
    """Adjust wave speed and time step to solve compatibility equations.
    
    Parameters
    ----------
    wn : wntr.network.model.WaterNetworkModel
        Network
    N : numpy array
        Number of discretization for each pipe
    
    Returns
    -------
    wn : wntr.network.model.WaterNetworkModel
        Network with adjusted wave speed.
    dt : float 
        Adjusted time step
    """

    from numpy import transpose as trans  

    phi = [np.longdouble(pipe.length / pipe.wavev / N[int(pipe.id)-1])
                            for _, pipe in wn.pipes()]
    phi = np.array(phi).reshape((len(phi), 1))
    theta = np.longdouble(1/ np.matmul(trans(phi), phi) * \
        np.matmul(trans(phi), np.ones((len(phi), 1))))
    
    # adjust time step 
    dt = np.float64(1/theta)

    # adjust the wave speed of each links    
    for _, pipe in wn.pipes():    
        pipe.wavev = np.float64(pipe.wavev * phi[int(pipe.id)-1] * theta)
    
    return wn, dt 

def discretization( wn, npipe, dt):
    """Discritize in temporal and spatial space using wave speed adjustement scheme. 
    
    Parameters
    ----------
    wn : wntr.network.model.WaterNetworkModel
        Network
    npipe : integer 
        Number of pipes
    dt : float 
        User defined time step

    Returns
    -------
    wn : wntr.network.model.WaterNetworkModel
        Network with updated parameters
    dt : float 
        Adjusted time step
    Ndis : numpy array 
        Number of discritization for each pipe
    """
    
    max_dt = max_time_step(wn)
    if dt > max_dt: 
        print("""time step is too large. Please define a time step that is
              less than %s """ %max_dt)
    else :
        Ndis = cal_N(wn, npipe, dt) 
        wn, dt = adjust_wavev(wn, Ndis)       
    return  wn, dt, Ndis

def valvesetting(dt, tf, valve_op):
    """Define valve operation curve (percentage open v.s. time)
    
    Parameters
    ----------
    dt : float
        Time step
    tf : float 
        Simulation Time
    valve_op : list 
        Contains paramtes to defie valve operation rule
        valve_op = [tc,ts,se,m]
        tc : the duration takes to close the valve [s]
        ts : closure start time [s]
        se : final open percentage [s]
        m  : closure constant [unitless]

    Returns
    -------
    s : list 
        valve operation curve 
    """

    [tc,ts,se,m] = valve_op
    tn = int(tf/dt)
    # aburupt closure 
    if tc ==0:
        s =  np.array([(1- (i*dt- ts))**1    for i in range(tn)])
        s[s>1] = 1
        s[s<1] = se
    # gradual closure 
    else:
        s =  np.array([(1- (i*dt- ts)/tc)**m    for i in range(tn)])
        s[s>1] = 1
        s[s<se] = se
    return s   

def pumpsetting(dt, tf, pump_op):
    """Define pump operation curve (percentage open v.s. time)
    
    Parameters
    ----------
    dt : float
        Time step
    tf : float 
        Simulation Time
    valve_op : list 
        Contains paramtes to defie valve operation rule
        valve_op = [tc,ts,se,m]
        tc : the duration takes to close the valve [s]
        ts : closure start time [s]
        se : final open percentage [s]
        m  : closure constant [unitless]

    Returns
    -------
    s : list 
        valve operation curve 
    """
    [tc,ts,se,m] = pump_op
    tn = int(tf/dt)
    # gradual closure 
    if tc != 0:
        s =  np.array([(1- (i*dt- ts)/tc)**m    for i in range(tn)])
        s[s>1] = 1
        s[s<se] = se
    
    # aburupt closure 
    if tc ==0:
        s =  np.array([(1- (i*dt- ts))**1    for i in range(tn)])
        s[s>1] = 1
        s[s<1] = se  
        
    return s