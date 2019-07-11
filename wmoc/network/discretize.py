"""
The wmoc.network.discretize contains methods to perform 
spatial and temporal discritization by adjusting wave speed 
and time step to solve compatibility equations in case of 
uneven wave travel time. 

"""

from __future__ import print_function
import numpy as np

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
    max_dt = np.inf
 
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
    N : numpy.ndarray
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
    Ndis : numpy.ndarray 
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
