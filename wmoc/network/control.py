"""
The wmoc.network.control module includes method to define
network controls of the pump and valve.These control modify
parameters in the network during trasient simulation.

"""
import numpy as np

def valveclosing(dt, tf, valve_op):
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

def valveopening(dt, tf, valve_op):
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
    # aburupt opening
    if tc ==0:
        s =  np.array([((i*dt- ts))**1    for i in range(tn)])
        s[s>0] = se
        s[s<0] = 0
    # gradual opening
    else:
        s =  np.array([((i*dt- ts)/tc)**m    for i in range(tn)])
        s[s<0] = 0
        s[s>se] = se
    return s

def pumpclosing(dt, tf, pump_op):
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
    # do not allow the pump to be fully cloased due to numerical issues
    if se == 0.:
        se = 0.0001

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

def pumpopening(dt, tf, pump_op):
    """Define pump operation curve (percentage open v.s. time)

    Parameters
    ----------
    dt : float
        Time step
    tf : float
        Simulation Time
    pump_op : list
        Contains paramtes to defie pump operation rule
        pump_op = [tc,ts,se,m]
        tc : the duration takes to start up the pump [s]
        ts : open start time [s]
        se : final open percentage [s]
        m  : closure constant [unitless]

    Returns
    -------
    s : list
        valve operation curve
    """

    [tc,ts,se,m] = pump_op
    tn = int(tf/dt)
    # aburupt opening
    if tc ==0:
        s =  np.array([((i*dt- ts))**1    for i in range(tn)])
        s[s>0] = se
        s[s<0] = 0
    # gradual opening
    else:
        s =  np.array([((i*dt- ts)/tc)**m    for i in range(tn)])
        s[s<0] = 0
        s[s>se] = se
    return s

def burstsetting(dt,tf,ts,tc,final_burst_coeff):
    """ Calculate the burst area as a function of simulation time

    Parameters
    ----------
    dt : float
        Time step
    tf : float
        Simulation Time
    ts : float
        Burst start time
    tc : float
        Time for burst to fully develop
    final_burst_coeff : list or float
        Final emitter coefficient at the burst nodes
    """

    tn = int(tf/dt)
    if tc !=0 :
        s = np.array([(i*dt- ts)/tc    for i in range(tn)])
        s[s>1] = 1
        s[s<0] = 0
        burst_A = final_burst_coeff * s
    else:
        s = np.array([(i*dt- ts)/tc    for i in range(tn)])
        s[s>1] = 1
        s[s<0] = 0
        burst_A = final_burst_coeff * s
    return burst_A
