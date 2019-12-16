"""
The tsnet.network.discretize contains methods to perform
spatial and temporal discretization by adjusting wave speed
and time step to solve compatibility equations in case of
uneven wave travel time.

"""

from __future__ import print_function
import numpy as np

def discretization(tm, dt):
    """Discretize in temporal and spatial space using wave speed adjustment scheme.

    Parameters
    ----------
    tm : tsnet.network.geometry.TransientModel
        Network
    dt : float
        User defined time step

    Returns
    -------
    tm : tsnet.network.geometry.TransientModel
        Network with updated parameters
    """

    max_dt = max_time_step(tm)
    if dt > max_dt:
        raise ValueError("time step is too large. Please define \
                    a time step that is less than %.5f " %max_dt)
    else :
        Ndis = cal_N(tm, dt)

        # add number of segments as a new attribute to each pipe
        i = 0
        for _, pipe in tm.pipes():
            pipe.number_of_segments = int(Ndis[i])
            i+=1
        # adjust wave speed and calculate time step
        tm = adjust_wavev(tm)
    return  tm


def max_time_step(tm):
    """Determine the maximum time step based on Courant's criteria.

    Parameters
    ----------
    tm : tsnet.network.geometry.TransientModel
        Network

    Returns
    -------
    max_dt : float
        Maximum time step allowed for this network
    """
    max_dt = np.inf

    for _, pipe in tm.pipes():
        dt = pipe.length / (2. * pipe.wavev)
        if max_dt > dt :
            max_dt = dt #- 0.001  # avoid numerical issue which cause N = 0
    return max_dt

def discretization_N(tm, dt):
    """Discretize in temporal and spatial space using wave speed adjustment scheme.

    Parameters
    ----------
    tm : tsnet.network.geometry.TransientModel
        Network
    dt : float
        User defined time step

    Returns
    -------
    tm : tsnet.network.geometry.TransientModel
        Network with updated parameters
    """

    Ndis = cal_N(tm, dt)

    # add number of segments as a new attribute to each pipe
    i = 0
    for _, pipe in tm.pipes():
        pipe.number_of_segments = int(Ndis[i])
        i+=1
    # adjust wave speed and calculate time step
    tm = adjust_wavev(tm)
    return  tm


def max_time_step_N(tm, N):
    """Determine the maximum time step based on Courant's criteria.

    Parameters
    ----------
    tm : tsnet.network.geometry.TransientModel
        Network

    Returns
    -------
    max_dt : float
        Maximum time step allowed for this network
    """
    max_dt = np.inf
    for _, pipe in tm.pipes():
        dt = pipe.length / (N * pipe.wavev)
        if max_dt > dt :
            max_dt = dt #- 1e-5  # avoid numerical issue which cause N = 0
    return max_dt

def cal_N(tm,  dt):
    """Determine the number of computation unites ($N_i$) for each pipes.

    Parameters
    ----------
    tm : tsnet.network.geometry.TransientModel
        Network
    dt : float
        Time step for transient simulation
    """
    N = np.zeros((tm.num_pipes,1))

    for _, pipe in tm.pipes() :
        # N[int(pipe.id)-1] =  int(2*np.int(pipe.length/ (2. * pipe.wavev *dt)))
        N[int(pipe.id)-1] =  round(np.int(pipe.length/ (pipe.wavev *dt)))
    return N


def adjust_wavev(tm):
    """Adjust wave speed and time step to solve compatibility equations.

    Parameters
    ----------
    tm : tsnet.network.geometry.TransientModel
        Network

    Returns
    -------
    tm : tsnet.network.geometry.TransientModel
        Network with adjusted wave speed.
    dt : float
        Adjusted time step
    """

    from numpy import transpose as trans
    phi = [np.longdouble(pipe.length / pipe.wavev / pipe.number_of_segments)
                            for _, pipe in tm.pipes()]
    phi = np.array(phi).reshape((len(phi), 1))
    tm.wavespeed_adj = np.sum(phi**2)
    theta = np.longdouble(1/ np.matmul(trans(phi), phi) * \
        np.matmul(trans(phi), np.ones((len(phi), 1))))

    # adjust time step
    dt = np.float64(1/theta)

    # adjust the wave speed of each links
    for _, pipe in tm.pipes():
        pipe.wavev = np.float64(pipe.wavev * phi[int(pipe.id)-1] * theta)

    # set time step as a new attribute to TransientModel
    tm.time_step =dt
    return tm

