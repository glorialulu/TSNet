"""
The wmoc.network.geometry read in the geometry defined by EPANet
.inp file, and assign aditional parameters needed in transient 
simution later in wmoc.

"""

from __future__ import print_function
import wntr
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import logging

logger = logging.getLogger(__name__)

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