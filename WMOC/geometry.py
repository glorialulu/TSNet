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
    """Take .inp fileas input and output the geometry of the network
    and parameters of the nodes and links. nodes and links are objects."""    
       
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
    for pn, pipe in wn.pipes():  
        pipe.wavev = lambda: None        
        setattr(pipe, 'wavev', wavev[i])
        i+=1
        
    # calculate the slope for each pipe 
    for pn, pipe in wn.pipes():  
        pipe.theta = lambda: None
        try: 
            theta = np.sin(np.arctan(pipe.end_node.elevation-pipe.start_node.elevation)/pipe.length)
        except:
            theta = 0.0
        setattr(pipe, 'theta',theta )
        
    # assign ID to each links, start from 1.      
    i =1
    for ln, link in wn.links():  
        link.id = lambda: None 
        setattr(link, 'id', i)
        i+=1 
   
    # assign ID to each links, start from 1.      
    i =1
    for nn, node in wn.nodes():  
        node.id = lambda: None 
        setattr(node, 'id', i)
        i+=1     ## Graph the network
        
    if plot: 
#        wntr.graphics.plot_network(wn, title=wn.name)
        # Compute betweenness centrality
        plt.figure(figsize=(4,4), dpi=80, facecolor='w', edgecolor='k') 
        bet_cen = nx.betweenness_centrality(G)
        wntr.graphics.plot_network(wn, node_attribute=bet_cen, node_size=30, 
                                title='Betweenness Centrality')
    return  wn, npipe

def topology(wn, npipe, pressure_zone_bc):
    G = wn.get_graph()
    
    length = wn.query_link_attribute('length')
    G.weight_graph(link_attribute = length)

    
    # add 'id' attribute to networkx links 
    i =1
    for ln, link in wn.links():          
        G.edges[link.start_node_name,link.end_node_name, ln]['id'] = i
        i+=1 
        
 
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
#        print ('Calculating pipe #',pipe.name)

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


        if  utype[pn][0] == 'Junction' and utype[pn][1].name in pressure_zone_bc:
            utype[pn] = ('PressureZone',utype[pn][1])
            
        if  dtype[pn][0] == 'Junction' and dtype[pn][1].name in pressure_zone_bc:
            dtype[pn] = ('PressureZone',dtype[pn][1])            

#        print('links1', links1[pn], 'links2', links2[pn])
#        print ('utype', utype[pn], 'dtype', dtype[pn])
        
    return links1, links2, utype, dtype

def max_time_step(wn):
    """Detrtmine the common maximum time step value for all pipes
    that satisfies Courant's criteria. 
    ."""
    max_dt = math.inf
 
    for name, pipe in wn.pipes():
        dt = pipe.length / (2. * pipe.wavev)
        if max_dt > dt :
            max_dt = dt           
    return max_dt

def cal_N(wn, npipe, dt):
    """ Determine the number of computation uites ($N_i$) for each pipes.""" 
    
    N = np.zeros((npipe,1))
    
    for _, pipe in wn.pipes() :
#        print ('N', int(pipe.id)-1)
        N[int(pipe.id)-1] =  int(2*np.int(pipe.length/ (2. * pipe.wavev *dt)))
    return N


def adjust_wavev( wn, N):
    """Adjust wave sppeed to solve the compatibility equations 
    in caseof an uneven wave travel time. 
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
#        adj = (phi[int(pipe.id)-1] * theta -1)*100
#        print ('wave speed of pipe %s ajusted  by %s ' %(pipe.id, adj), '%')
        pipe.wavev = np.float64(pipe.wavev * phi[int(pipe.id)-1] * theta)
    
    return wn, dt 

def discretization( wn, npipe, dt ):
    """ Discritize in temporal and spatial space. 
    Caculate the time interval (dt) using wave speed adjustement scheme 
    to solve the compatibility equations in caseof an uneven wave travel time."""
    
    max_dt = max_time_step(wn)
    if dt > max_dt: 
        print("""time step is too large. Please define a time step that is
              less than %s """ %max_dt)
    else :
        Ndis = cal_N(wn, npipe, dt) 
        wn, dt = adjust_wavev(wn, Ndis)

#        for _, pipe in wn.pipes() :
#            print ('pipe', pipe.id, 'spatial discretization', Ndis[int(pipe.id)-1],
#                    'wave spped', pipe.wavev, 
#                    'dt', pipe.length/Ndis[int(pipe.id)-1]/ pipe.wavev)
        
    return  wn, dt, Ndis

def valvesetting(dt, tf, valve_op):
    "Define valve operation curve."
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
    """ Define pump operation curve. """
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
