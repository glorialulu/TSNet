"""
The wmoc.simulation.solver module contains methods to solver MOC
for different grid configurations, including:
1. inner_node
2. valve_node
3. pump_node
4. source_pump
5. valve_end
6. dead_end
7. rev_end
8. add_leakage

"""
from __future__ import print_function
import numpy as np
import warnings

def inner_node(link1, link2, demand, H1, V1, H2, V2, dt, g, nn, s1, s2):
    """Inner boudary MOC using C+ and C- characteristic curve

    Parameters
    ----------
    link1 : object
        Pipe object of C+ charateristics curve
    link2 : object
        Pipe object of C- charateristics curve
    demand : float
        demand at the junction
    H1 : list
        List of the head of C+ charateristics curve
    V1 : list
        List of the velocity of C+ charateristics curve
    H2 : list
        List of the head of C- charateristics curve
    V2 : list
        List of the velocity of C- charateristics curve
    dt : float
        Time step
    g : float
        Gravity aceleration
    nn : int
        The index of the calculation node
    s1 : list
        List of signs that represent the direction of the flow
        in C+ charateristics curve
    s2 : list
        List of signs that represent the direction of the flow
        in C- charateristics curve

    Returns
    -------
    HP : float
        Head at current node at current time
    VP : float
        Velocity at current node at current time
    """

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
    A1 = [np.pi * D1[i]**2. / 4.  for i in range(len(link1))]   # m^2
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
    A2 = [np.pi * D2[i]**2. / 4.  for i in range(len(link2))]    #m^2
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


def valve_node(KL_inv, link1, link2, H1, V1, H2, V2, dt, g, nn, s1, s2):
    """Inline valve node MOC calculation

    Parameters
    ----------
    KL_inv : int
        Inverse of the valve loss coefficient at current time
    link1 : object
        Pipe object of C+ charateristics curve
    link2 : object
        Pipe object of C- charateristics curve
    H1 : list
        List of the head of C+ charateristics curve
    V1 : list
        List of the velocity of C+ charateristics curve
    H2 : list
        List of the head of C- charateristics curve
    V2 : list
        List of the velocity of C- charateristics curve
    dt : float
        Time step
    g : float
        Gravity aceleration
    nn : int
        The index of the calculation node
    s1 : list
        List of signs that represent the direction of the flow
        in C+ charateristics curve
    s2 : list
        List of signs that represent the direction of the flow
        in C- charateristics curve
    """

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
    A1 = [np.pi * D1[i]**2. / 4.  for i in range(len(link1))]   #m^2
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
    A2 = [np.pi * D2[i]**2. / 4.  for i in range(len(link2))]    #m^2
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
        warnings.warn('Error: The quadratic equation has no real solution (valve)')

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
            warnings.warn('Error: The quadratic equation has no real solution (valve)')

        if nn == 0:  # pipe start
            VP = VP*A1[0]/A2[0]
            HP = (C2[0,0] + VP ) / C2[0,1]
        else:        # pipe end
            VP = VP
            HP = (C1[0,0] - VP) / C1[0,1]
    return HP, VP


def pump_node(pumpc,link1, link2, H1, V1, H2, V2, dt, g, nn, s1, s2):
    """ Inline pump node MOC calculation

    Parameters
    ----------
    pumpc : list
        Parameters (a, b,c) to define pump characteristic cure,
        so that
        .. math:: h_p = a*Q**2 + b*Q + c
    link1 : object
        Pipe object of C+ charateristics curve
    link2 : object
        Pipe object of C- charateristics curve
    H1 : list
        List of the head of C+ charateristics curve
    V1 : list
        List of the velocity of C+ charateristics curve
    H2 : list
        List of the head of C- charateristics curve
    V2 : list
        List of the velocity of C- charateristics curve
    dt : float
        Time step
    g : float
        Gravity acceleration
    nn : int
        The index of the calculation node
    s1 : list
        List of signs that represent the direction of the flow
        in C+ charateristics curve
    s2 : list
        List of signs that represent the direction of the flow
        in C- charateristics curve
    """

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
    A1 = [np.pi * D1[i]**2. / 4.  for i in range(len(link1))] # m^2
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
    A2 = [np.pi * D2[i]**2. / 4.  for i in range(len(link2))]# m^2
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
        warnings.warn('Error: The quadratic equation has no real solution (pump)')

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
    """Source Pump boundary MOC calculation

    Parameters
    ----------
    pump : list
        pump[0]: elevation of the reservoir/tank
        pump[1]: Parameters (a, b,c) to define pump characteristic cure,
        so that
        .. math:: h_p = a*Q**2 + b*Q + c
    link2 : object
        Pipe object of C- charateristics curve
    H2 : list
        List of the head of C- charateristics curve
    V2 : list
        List of the velocity of C- charateristics curve
    dt : float
        Time step
    g : float
        Gravity acceleration
    s2 : list
        List of signs that represent the direction of the flow
        in C- charateristics curve
    """
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
    A2 = [np.pi * D2[i]**2. / 4.  for i in range(len(link2))]    #m^2
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

    # parameters of the quadratic polynomial
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
        warnings.warn('The quadratic equation has no real solution (pump)')

    if HP > Hsump:
        VP = np.float64(-C2[0,0] + C2[0,1]*HP)
    else :
        HP = Hsump
        VP = np.float64(-C2[0,0] + C2[0,1]*HP)

    if VP <= 0 : # positive flow
        warnings.warn( "Reverse flow stopped by check valve!")
        VP = 0
        HP = (C2[0,0] + VP ) / C2[0,1]

    return HP, VP



def valve_end(H1, V1, V, nn, a, g, f, D, dt):
    """ End Valve boundary MOC calculation

    Parameters
    ----------
    H1 : float
        Head of the C+ charateristics curve
    V1 : float
        Velocity of the C+ charateristics curve
    V : float
        Velocity at the valve end at current time
    nn : int
        The index of the calculation node
    a : float
        Wave speed at the valve end
    g : float
        Gravity acceleration
    f : float
        friction factor of the current pipe
    D : float
        diameter of the current pipe
    dt : float
        Time step
    """
    if nn == 0 :
        HP = H1 + a/g*(V - V1) + a/g*f*dt/(2.*D)*V1*abs(V1)
        VP = V
    else :
        HP = H1 - a/g*(V - V1) - a/g*f*dt/(2.*D)*V1*abs(V1)
        VP = V
    return HP,VP

def dead_end(linkp, H1, V1, nn, a, g, f, D, dt):
    """Dead end boundary MOC calculation with pressure dependant demand

    Parameters
    ----------
    link1 : object
        Current pipe
    H1 : float
        Head of the C+ charateristics curve
    V1 : float
        Velocity of the C+ charateristics curve
    nn : int
        The index of the calculation node
    a : float
        Wave speed at the valve end
    g : float
        Gravity acceleration
    f : float
        friction factor of the current pipe
    D : float
        diameter of the current pipe
    dt : float
        Time step
    """

    A = np.pi/4. * linkp.diameter**2.

    if nn == 0:
        k = linkp.start_node.demand_coeff
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
            warnings.warn("The quadratic equation has no real solution (dead end).\
The resuls might not be accurate.")
        VP = V1 - g/a *H1 - f*dt/(2.*D)*V1*abs(V1) + g/a*HP
    else :
        k = linkp.end_node.demand_coeff
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
            warnings.warn("The quadratic equation has no real solution (dead end).\
The resuls might not be accurate.")
        VP = V1 + g/a *H1 - f*dt/(2.*D)*V1*abs(V1) - g/a*HP
    return HP,VP

def rev_end( H2, V2, H, nn, a, g, f, D, dt):
    """Reservoir/ Tank boundary MOC calculation

    Parameters
    ----------
    H2 : list
        List of the head of C- charateristics curve
    V2 : list
        List of the velocity of C- charateristics curve
    H : float
        Head of the reservoir/tank
    nn : int
        The index of the calculation node
    a : float
        Wave speed at the valve end
    g : float
        Gravity acceleration
    f : float
        friction factor of the current pipe
    D : float
        diameter of the current pipe
    dt : float
        Time step
    """
    if nn == 0 :
        VP = V2 + g/a*(H - H2) - f*dt/(2.*D)*V2*abs(V2)
        HP = H
    else:
        VP = V2 - g/a*(H - H2) - f*dt/(2.*D)*V2*abs(V2)
        HP = H
    return HP, VP

def add_leakage(emitter_coef, link1, link2, H1, V1, H2, V2, dt, g, nn, s1, s2):
    r"""Leakge Node MOC calculation

    Parameters
    ----------
    emitter_coef : flot
        float, optional
        Required if leak_loc is defined
        The leakage coefficient of the leakge
        .. math:: Q_leak = leak_A  [ m^3/s/(m H20)^(1/2)] * \sqrt(H)
    link1 : object
        Pipe object of C+ charateristics curve
    link2 : object
        Pipe object of C- charateristics curve
    H1 : list
        List of the head of C+ charateristics curve
    V1 : list
        List of the velocity of C+ charateristics curve
    H2 : list
        List of the head of C- charateristics curve
    V2 : list
        List of the velocity of C- charateristics curve
    dt : float
        Time step
    g : float
        Gravity aceleration
    nn : int
        The index of the calculation node
    s1 : list
        List of signs that represent the direction of the flow
        in C+ charateristics curve
    s2 : list
        List of signs that represent the direction of the flow
        in C- charateristics curve
    """

    emitter_coef = emitter_coef  # m^3/s//(m H2O)^(1/2)

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
    A1 = [np.pi * D1[i]**2. / 4.  for i in range(len(link1))]   #m^2
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
    A2 = [np.pi * D2[i]**2. / 4.  for i in range(len(link2))]    #m^2
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
        warnings.warn('Error: The quadratic equation has no real solution (leakage)')

    if nn == 0:  # pipe start
        VP = np.float64(-C2[:,0]+ C2[:,1]*HP)
    else:        # pipe end
        VP = np.float64(C1[:,0] - C1[:,1]*HP)
    return HP, VP