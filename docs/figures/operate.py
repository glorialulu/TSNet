import numpy as np
import matplotlib.pyplot as plt
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

dt = 0.1
tf = 2
t = np.arange(0, tf, dt)

tc = 1 # pump closure period
ts = 0 # pump closure start time
se = 0 # end open percentage

m = 1 # closure constant
valve_op = [tc,ts,se,m]
s1 = valveclosing(dt, tf, valve_op)

m = 2 # closure constant
valve_op = [tc,ts,se,m]
s2 = valveclosing(dt, tf, valve_op)

tc = 1 # pump closure period
ts = 0 # pump closure start time
se = 1 # end open percentage

m = 1 # closure constant
valve_op = [tc,ts,se,m]
s3 = valveopening(dt, tf, valve_op)

m = 2 # closure constant
valve_op = [tc,ts,se,m]
s4 = valveopening(dt, tf, valve_op)

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
fig = plt.figure(figsize=(10,4), dpi=80, facecolor='w', edgecolor='k')
plt.subplot(1, 2, 1)
plt.plot(t,s1,'k-', label='m=1')
plt.plot(t,s2,'r-', label='m=2')
plt.xlim([t[0],t[-1]])
plt.title('(a)')
plt.xticks([0,1,t[-1]],['ts','ts + tc','tf'])
plt.yticks([0,0.2,0.4,0.6,0.8,1.0],['se','0.2','0.4','0.6','0.8','1.0'])
plt.xlabel("Time [s]")
plt.ylabel("Valve opening ratio")
plt.legend(loc='best')
plt.grid(True)

plt.subplot(1, 2, 2)
plt.plot(t,s3,'k-', label='m=1')
plt.plot(t,s4,'r-', label='m=2')
plt.xlim([t[0],t[-1]])
plt.title('(b)')
plt.xticks([0,1,t[-1]],['ts','ts + tc','tf'])
plt.yticks([0,0.2,0.4,0.6,0.8,1.0],['0','0.2','0.4','0.6','0.8','se'])
plt.xlabel("Time [s]")
plt.legend(loc='best')
plt.grid(True)

fig.savefig('./docs/figures/valve_operating.pdf', format='pdf',dpi=500)
