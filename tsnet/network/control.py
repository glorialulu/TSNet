"""
The tsnet.network.control module includes method to define
network controls of the pump and valve.These control modify
parameters in the network during transient simulation.

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
        Contains parameters to define valve operation rule
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
    # abrupt closure
    if tc ==0:
        s =  np.array([(1- (i*dt- ts))**1    for i in range(tn)])
        s[s>1] = 1
        s[s<1] = se
    # gradual closure
    else:
        t = np.array([(i*dt- ts)/tc for i in range(tn)])
        t[t>1] = 1
        t[t<0] = 0
        s =  np.array([1 - (1-se)*t[i]**m for i in range(tn)])
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
        Contains parameters to define valve operation rule
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
    # abrupt opening
    if tc ==0:
        s =  np.array([((i*dt- ts))**1    for i in range(tn)])
        s[s>0] = se
        s[s<0] = 0
    # gradual opening
    else:
        t = np.array([(i*dt- ts)/tc for i in range(tn)])
        t[t>1] = 1
        t[t<0] = 0
        s =  np.array([se* (t[i])**m for i in range(tn)])
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
        Contains parameters to define valve operation rule
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
    # do not allow the pump to be fully closed due to numerical issues
    if se == 0.:
        se = 0.0001

    tn = int(tf/dt)
    # gradual closure
    if tc != 0:
        s =  np.array([(1- (i*dt- ts)/tc)**m    for i in range(tn)])
        s[s>1] = 1
        s[s<se] = se

    # abrupt closure
    if tc ==0:
        t = np.array([(i*dt- ts)/tc for i in range(tn)])
        t[t>1] = 1
        t[t<0] = 0
        s =  np.array([1 - (1-se)*t[i]**m for i in range(tn)])
        s[s>1] = 1
        s[s<se] = se
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
        Contains parameters to define pump operation rule
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
    # abrupt opening
    if tc ==0:
        s =  np.array([((i*dt- ts))**1    for i in range(tn)])
        s[s>0] = se
        s[s<0] = 0
    # gradual opening
    else:
        t = np.array([(i*dt- ts)/tc for i in range(tn)])
        t[t>1] = 1
        t[t<0] = 0
        s =  np.array([se* (t[i])**m for i in range(tn)])
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
        s = np.array([(i*dt- ts)    for i in range(tn)])
        s[s>0] = 1
        s[s<0] = 0
        burst_A = final_burst_coeff * s

    return burst_A


def demandpulse(dt, tf, tc, ts, tp, dp):
    """Calculate demand pulse multiplier

    Parameters
    dt : float
        Time step
    tf : float
        Simulation Time
    tc : float
        Total pulse duration
    ts : float
        Pulse start time
    tp : float
        Pulse increase time
    dp : float
        Pulse multiplier
    """
    tn = int(tf/dt)
    x = np.linspace(0,tf,tn)
    t_change = tp
    stay = tc - 2*tp
    if t_change !=0 :
        s = np.piecewise(x, [x<=ts, (x>ts)& (x<= ts+t_change),
                        (x>ts+t_change) & (x<=ts+t_change+stay),
                        (x>ts+t_change+stay) & (x<=ts+tc),
                        x>ts+tc],
                        [0, lambda x: (x-ts)/t_change,
                        1, lambda x: 1- (x-ts-t_change-stay)/t_change,
                        0])

        pulse_mul = dp * s

    else:
        s = np.piecewise(x, [x<=ts,(x>ts) & (x<=ts+stay), x>ts+tc],
                        [0, 1, 0])
        pulse_mul = dp * s

    # import matplotlib.pyplot as plt

    # fig1 = plt.figure(figsize=(4,4), dpi=150, facecolor='w', edgecolor='k')
    # plt.plot(x,pulse_mul, 'k', lw=2.5)
    # plt.xlim(ts-0.2,ts+tc+0.2)
    # plt.ylim(-0.1,1.2)
    # plt.xticks([ts,ts+tp,ts+tc],('ts', 'ts+tp', 'ts+tc'))
    # plt.yticks([dp], ['dp'])
    # plt.xlabel("Time [s]")
    # plt.ylabel("Demand pulse multiplier")
    # plt.vlines(ts+tp, -0.1, 1, 'k', linestyles='dotted')
    # plt.vlines(ts, -0.1, 0, 'k', linestyles='dotted')
    # plt.vlines(ts+tc, -0.1, 0, 'k', linestyles='dotted')
    # plt.hlines(dp, ts-0.2, ts+tp, 'k', linestyles='dotted')
    # plt.show()
    # fig1.savefig("DemandMultiplier.pdf", format='pdf',dpi=500, bbox_inches = 'tight')
    return pulse_mul
