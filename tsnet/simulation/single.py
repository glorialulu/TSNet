"""
The tsnet.simulation.single contains methods to perform MOC
transient simulation on a single pipe, including
1. inner pipe
2. left boundary pipe (without C- charateristic grid)
3. right boundary pipe (without C+ characteristic grid)

"""
import numpy as np
from tsnet.simulation.solver import (
    inner_node_steady,
    inner_node_quasisteady,
    inner_node_unsteady,
    valve_node,
    pump_node,
    source_pump,
    valve_end,
    dead_end,
    rev_end,
    add_leakage,
    surge_tank,
    air_chamber
)

def inner_pipe (linkp, pn, dt, links1, links2, utype, dtype, p,
                H0, V0, H, V, H10, V10, H20, V20, pump, valve,
                friction, dVdt, dVdx,
                dVdt10, dVdx10, dVdt20, dVdx20):
    """MOC solution for an individual inner pipe.

    Parameters
    ----------
    linkp : object
        Current pipe object
    pn : int
        Current pipe ID
    dt : float
        Time step
    H : numpy.ndarray
        Head of current pipe at current time step [m]
    V : numpy.ndarray
        Velocity of current pipe at current time step [m/s]
    links1 : list
        Upstream adjacent pipes
    links2 : list
        Downstream adjacent pipes
    utype : list
        Upstream adjacent link type, and if not pipe, their name
    dtype : list
        Downstream adjacent link type, and if not pipe, their name
    p : list
        pipe list
    H0 : numpy.ndarray
        Head of current pipe at previous time step [m]
    V0 : numpy.ndarray
        Velocity of current pipe at previous time step [m/s]
    H10 : list
        Head of left adjacent nodes at previous time step [m]
    V10 : list
        Velocity of left adjacent nodes at previous time step [m/s]
    H20 : list
        Head of right adjacent nodes at previous time step [m]
    V20 : list
        Velocity of right adjacent nodes at previous time step [m/s]
    pump : list
        Characteristics of the pump
    valve : list
        Characteristics of the valve
    friction: str
        friction model, e.g., 'steady', 'quasi-steady', 'unsteady',
        by default 'steady'
    dVdt: numpy.ndarray
        local instantaneous acceleration approximation to be used
        for unsteady friction calculation, 0 if not
        in unsteady friction mode [m/s^2]
    dVdx: numpy.ndarray
        convective instantaneous acceleration approximation to be used
        for unsteady friction calculation, 0 if not
        in unsteady friction mode [m/s^2]
    dVdt10 : list
        local instantaneous acceleration of left adjacent nodes at previous time step [m]
    dVdx10 : list
        convective instantaneous acceleration of left adjacent nodes at previous time step [m/s]
    dVdt20 : list
        local instantaneous acceleration of right adjacent nodes at previous time step [m]
    dVdx20 : list
        convective instantaneous acceleration of right adjacent nodes at previous time step [m/s]

    Returns
    -------
    H : numpy.ndarray
        Head results of the current pipe at current time step. [m]
    V : numpy.ndarray
        Velocity results of the current pipe at current time step. [m/s]
    """

    # Properties of current pipe
    g = 9.8                          # m/s^2
    link1 = [p[abs(i)-1] for i in links1]
    link2 = [p[abs(i)-1] for i in links2]
    n = linkp.number_of_segments    # spatial discretization

    # inner nodes
    if friction == 'steady':
        H[1:-1], V[1:-1] = inner_node_steady(linkp, H0, V0, dt, g)
    elif friction == 'quasi-steady':
        H[1:-1], V[1:-1] = inner_node_quasisteady(linkp, H0, V0, dt, g)
    else:
        H[1:-1], V[1:-1] = inner_node_unsteady(linkp, H0, V0, dt, g,
             dVdx, dVdt)

    # Pipe start
    V1 = V10;     H1 = H10       #list
    V2 = V0[1];   H2 = H0[1]
    dVdx1 = dVdx10 ; dVdt1 = dVdt10
    dVdx2 =  dVdx[0]; dVdt2 = dVdt[1]

    if utype[0] == 'Pipe':
        if linkp.start_node.transient_node_type == 'SurgeTank':
            shape = linkp.start_node.tank_shape
            H[0], V[0], Qs = surge_tank(shape, link1, linkp,
                H1, V1, H2, V2, dt, g, 0,  np.sign(links1), [-1],
                friction, dVdx1, dVdx2, dVdt1, dVdt2)
            linkp.start_node.water_level = H[0]
            linkp.start_node.tank_flow = Qs
        elif linkp.start_node.transient_node_type == 'Chamber':
            shape = linkp.start_node.tank_shape
            H[0], V[0], Qs, zp = air_chamber(shape, link1, linkp,
                H1, V1, H2, V2, dt, g, 0,  np.sign(links1), [-1],
                friction, dVdx1, dVdx2, dVdt1, dVdt2)
            linkp.start_node.water_level = zp
            linkp.start_node.tank_flow = Qs
        else:
            elev = linkp.start_node.elevation
            emitter_coeff = linkp.start_node.emitter_coeff + linkp.start_node.demand_coeff
            block_per = linkp.start_node.block_per
            H[0], V[0] = add_leakage(emitter_coeff, block_per, link1, linkp, elev,
                H1, V1, H2, V2, dt, g, 0,  np.sign(links1), [-1],
                friction, dVdx1, dVdx2, dVdt1, dVdt2)
    elif utype[0] == 'Pump':
        pumpc = pump[0]
        H[0], V[0] = pump_node(pumpc, link1, linkp,
            H1, V1, H2, V2, dt, g, 0,  np.sign(links1), [-1],
            friction, dVdx1, dVdx2, dVdt1, dVdt2)
    elif utype[0] == 'Valve':
        valvec = valve[0]
        H[0], V[0] = valve_node(valvec, link1, linkp,
            H1, V1, H2, V2, dt, g, 0,  np.sign(links1), [-1],
            friction, dVdx1, dVdx2, dVdt1, dVdt2)

    # Pipe end
    V1 = V0[n-1];    H1 = H0[n-1]
    V2 = V20;        H2 = H20
    dVdx1 = dVdx[n-1] ; dVdt1 = dVdt[n-1]
    dVdx2 =  dVdx20; dVdt2 = dVdt20
    if dtype[0] == 'Pipe':
        if linkp.end_node.transient_node_type == 'SurgeTank':
            shape = linkp.end_node.tank_shape
            H[n], V[n], Qs = surge_tank(shape, linkp, link2,
                H1, V1, H2, V2, dt, g, n, [1], np.sign(links2),
                friction, dVdx1, dVdx2, dVdt1, dVdt2)
            linkp.end_node.water_level = H[n]
            linkp.end_node.tank_flow = Qs
        elif linkp.end_node.transient_node_type == 'Chamber':
            shape = linkp.end_node.tank_shape
            H[n], V[n], Qs,zp = air_chamber(shape, linkp, link2,
                H1, V1, H2, V2, dt, g, n, [1], np.sign(links2),
                friction, dVdx1, dVdx2, dVdt1, dVdt2)
            linkp.end_node.water_level = zp
            linkp.end_node.tank_flow = Qs
        else:
            elev = linkp.end_node.elevation
            emitter_coeff = linkp.end_node.emitter_coeff + linkp.end_node.demand_coeff
            block_per =  linkp.end_node.block_per
            H[n], V[n] = add_leakage(emitter_coeff, block_per,linkp, link2, elev,
                H1, V1, H2, V2, dt, g, n, [1], np.sign(links2),
                friction, dVdx1, dVdx2, dVdt1, dVdt2)
    elif dtype[0] == 'Pump':
        pumpc = pump[1]
        H[n], V[n] = pump_node(pumpc, linkp, link2,
                H1, V1, H2, V2, dt, g, n, [1], np.sign(links2),
                friction, dVdx1, dVdx2, dVdt1, dVdt2)

    elif dtype[0] == 'Valve':
        valvec = valve[1]
        H[n], V[n] = valve_node(valvec, linkp, link2,
                H1, V1, H2, V2, dt, g, n, [1], np.sign(links2),
                friction, dVdx1, dVdx2, dVdt1, dVdt2)
    return H, V

def left_boundary(linkp, pn, H, V, H0, V0, links2, p, pump, valve, dt,
                    H20, V20, utype, dtype,
                    friction, dVdt, dVdx, dVdt20, dVdx20) :
    """MOC solution for an individual left boundary pipe.

    Parameters
    ----------
    linkp : object
        Current pipe object
    pn : int
        Current pipe ID
    H : numpy.ndarray
        Head of current pipe at current time step [m]
    V : numpy.ndarray
        Velocity of current pipe at current time step [m/s]
    links2 : list
        Downstream adjacent pipes
    p : list
        pipe list
    pump : list
        Characteristics of the pump
    valve : list
        Characteristics of the valve
    n : int
        Number of discretization of current pipe
    dt : float
        Time step
    H0 : numpy.ndarray
        Head of current pipe at previous time step [m]
    V0 : numpy.ndarray
        Velocity of current pipe at previous time step [m/s]
    H20 : list
        Head of right adjacent nodes at previous time step [m]
    V20 : list
        Velocity of right adjacent nodes at previous time step [m/s]
    utype : list
        Upstream adjacent link type, and if not pipe, their name
    dtype : list
        Downstream adjacent link type, and if not pipe, their name
    friction: str
        friction model, e.g., 'steady', 'quasi-steady', 'unsteady',
        by default 'steady'
    dVdt: numpy.ndarray
        local instantaneous velocity approximation to be used
        for unsteady friction calculation, 0 if not
        in unsteady friction mode [m/s^2]
    dVdx: numpy.ndarray
        convective instantaneous velocity approximation to be used
        for unsteady friction calculation, 0 if not
        in unsteady friction mode [m/s^2]
    dVdt20 : list
        local instantaneous acceleration of right adjacent nodes at previous time step [m]
    dVdx20 : list
        convective instantaneous acceleration of right adjacent nodes at previous time step [m/s]

    Returns
    -------
    H : numpy.ndarray
        Head results of the current pipe at current time step. [m]
    V : numpy.ndarray
        Velocity results of the current pipe at current time step. [m/s]
    """

    link2 = [p[abs(i)-1] for i in links2]
    # Properties of current pipe
    f = linkp.roughness              # unitless
    D = linkp.diameter               # m
    g = 9.8                          # m/s^2
    a = linkp.wavev    # m/s
    n = linkp.number_of_segments   # spatial discretization
    KD = linkp.roughness_height

    # inner nodes
    if friction == 'steady':
        H[1:-1], V[1:-1] = inner_node_steady(linkp, H0, V0, dt, g)
    elif friction == 'quasi-steady':
        H[1:-1], V[1:-1] = inner_node_quasisteady(linkp, H0, V0, dt, g)
    else:
        H[1:-1], V[1:-1] = inner_node_unsteady(linkp, H0, V0, dt, g,
             dVdx, dVdt)

    # Pipe start (outer boundayr conditions)
    V2 = V0[1]; H2 = H0[1]
    dVdx2 = dVdx[0]; dVdt2= dVdt[1]
    if utype[0] == 'Reservoir' or  utype[0] == 'Tank':
        H[0], V[0] = rev_end (H2, V2, H[0], 0, a, g, f, D, dt,
        KD, friction, dVdx2, dVdt2)
    elif utype[0] == 'Valve':
        H[0], V[0] = valve_end (H2, V2, V[0], 0, a, g, f, D, dt,
        KD, friction, dVdx2, dVdt2)
    elif utype[0] == 'Junction':
        elev = linkp.start_node.elevation
        H[0], V[0] = dead_end (linkp , H2, V2, elev, 0, a, g, f, D, dt,
        KD, friction, dVdx2, dVdt2)
    elif utype[0] == 'Pump':  #source pump
        H[0], V[0] = source_pump(pump[0], linkp, H2, V2, dt, g, [-1],
        friction, dVdx2, dVdt2)

    # Pipe end  (inner boundary conditions)
    V1 = V0[n-1]; H1 = H0[n-1]     # upstream node
    V2 = V20;     H2 = H20         # downstream nodes
    dVdx1 = dVdx[n-1] ; dVdx2 = dVdx20
    dVdt1 = dVdt[n-1] ; dVdt2 = dVdt20

    if dtype[0] == 'Pipe':
        if linkp.end_node.transient_node_type == 'SurgeTank':
            shape = linkp.end_node.tank_shape
            H[n], V[n], Qs = surge_tank(shape, linkp, link2,
                H1, V1, H2, V2, dt, g, n, [1], np.sign(links2),
                friction, dVdx1, dVdx2, dVdt1, dVdt2)
            linkp.end_node.water_level = H[n]
            linkp.end_node.tank_flow = Qs

        elif linkp.end_node.transient_node_type == 'Chamber':
            shape = linkp.end_node.tank_shape
            H[n], V[n], Qs, zp = air_chamber(shape, linkp, link2,
                H1, V1, H2, V2, dt, g, n, [1], np.sign(links2),
                friction, dVdx1, dVdx2, dVdt1, dVdt2)
            linkp.end_node.water_level = zp
            linkp.end_node.tank_flow = Qs
        else:
            elev = linkp.end_node.elevation
            emitter_coeff = linkp.end_node.emitter_coeff + linkp.end_node.demand_coeff
            block_per =  linkp.end_node.block_per
            H[n], V[n] = add_leakage(emitter_coeff, block_per,linkp, link2, elev,
                H1, V1, H2, V2, dt, g, n, [1], np.sign(links2),
                friction, dVdx1, dVdx2, dVdt1, dVdt2)

    elif dtype[0] == 'Pump':
        pumpc = pump[1]
        H[n], V[n] = pump_node(pumpc, linkp, link2,
                H1, V1, H2, V2, dt, g, n, [1], np.sign(links2),
                friction, dVdx1, dVdx2, dVdt1, dVdt2)

    elif dtype[0] == 'Valve':
        valvec = valve[1]
        if links2 == []:
            H[i], V[i] = valve_end (H1, V1, V[i], n, a, g, f, D, dt,
            KD, friction, dVdx1, dVdt1)
        else:
            H[n], V[n] = valve_node(valvec, linkp, link2,
                H1, V1, H2, V2, dt, g, n, [1], np.sign(links2),
                friction, dVdx1, dVdx2, dVdt1, dVdt2)

    elif dtype[0] == 'Junction':
        elev = linkp.end_node.elevation
        H[n], V[n] = dead_end (linkp, H1, V1, elev, n, a, g, f, D, dt,
        KD, friction, dVdx1, dVdt1)

    return H, V

def right_boundary(linkp, pn, H0, V0, H, V, links1, p, pump, valve, dt,
                 H10, V10, utype, dtype,
                 friction, dVdt, dVdx, dVdt10, dVdx10):
    """MOC solution for an individual right boundary pipe.

    Parameters
    ----------
    linkp : object
        Current pipe object
    pn : int
        Current pipe ID
    H : numpy.ndarray
        Head of current pipe at current time step [m]
    V : numpy.ndarray
        Velocity of current pipe at current time step [m/s]
    links1 : list
        Upstream adjacent pipes
    p : list
        pipe list
    pump : list
        Characteristics of the pump
    valve : list
        Characteristics of the valve
    n : int
        Number of discretization of current pipe
    dt : float
        Time step
    H0 : numpy.ndarray
        Head of current pipe at previous time step [m]
    V0 : numpy.ndarray
        Velocity of current pipe at previous time step [m/s]
    H10 : list
        Head of left adjacent nodes at previous time step [m]
    V10 : list
        Velocity of left adjacent nodes at previous time step [m/s]
    utype : list
        Upstream adjacent link type, and if not pipe, their name
    dtype : list
        Downstream adjacent link type, and if not pipe, their name
    friction: str
        friction model, e.g., 'steady', 'quasi-steady', 'unsteady',
        by default 'steady'
    dVdt: numpy.ndarray
        local instantaneous velocity approximation to be used
        for unsteady friction calculation, 0 if not
        in unsteady friction mode [m/s^2]
    dVdx: numpy.ndarray
        convective instantaneous velocity approximation to be used
        for unsteady friction calculation, 0 if not
        in unsteady friction mode [m/s^2]
    dVdt10 : list
        local instantaneous acceleration of left adjacent nodes at previous time step [m]
    dVdx10 : list
        convective instantaneous acceleration of left adjacent nodes at previous time step [m/s]

   Returns
    -------
    H : numpy.ndarray
        Head results of the current pipe at current time step. [m]
    V : numpy.ndarray
        Velocity results of the current pipe at current time step. [m/s]
    """

    # Properties of current pipe
    link1 = [p[abs(i)-1] for i in links1]
    f = linkp.roughness              # unitless
    D = linkp.diameter               # m
    g = 9.8                          # m/s^2
    a = linkp.wavev                  # m/s
    n = linkp.number_of_segments   # spatial discretization
    KD = linkp.roughness_height

    # inner nodes
    if friction == 'steady':
        H[1:-1], V[1:-1] = inner_node_steady(linkp, H0, V0, dt, g)
    elif friction == 'quasi-steady':
        H[1:-1], V[1:-1] = inner_node_quasisteady(linkp, H0, V0, dt, g)
    else:
        H[1:-1], V[1:-1] = inner_node_unsteady(linkp, H0, V0, dt, g,
             dVdx, dVdt)

    # Pipe start (inner boundary conditions)
    V1 = V10; H1 = H10            # upstream node
    V2 = V0[1]; H2 = H0[1]    # downstream node
    dVdx1 = dVdx10 ; dVdx2 = dVdx[0]
    dVdt1 = dVdt10 ; dVdt2 = dVdt[1]
    if utype[0] == 'Pipe':
        if linkp.start_node.transient_node_type == 'SurgeTank':
            shape = linkp.start_node.tank_shape
            H[0], V[0], Qs = surge_tank(shape, link1, linkp,
                H1, V1, H2, V2, dt, g, 0,  np.sign(links1), [-1],
                friction, dVdx1, dVdx2, dVdt1, dVdt2)
            linkp.start_node.water_level = H[0]
            linkp.start_node.tank_flow = Qs
        if linkp.start_node.transient_node_type == 'Chamber':
            shape = linkp.start_node.tank_shape
            H[0], V[0], Qs, zp = air_chamber(shape, link1, linkp,
                H1, V1, H2, V2, dt, g, 0,  np.sign(links1), [-1],
                friction, dVdx1, dVdx2, dVdt1, dVdt2)
            linkp.start_node.water_level = zp
            linkp.start_node.tank_flow = Qs

        else:
            elev = linkp.start_node.elevation
            emitter_coeff = linkp.start_node.emitter_coeff + linkp.start_node.demand_coeff
            block_per =  linkp.start_node.block_per
            H[0], V[0] = add_leakage(emitter_coeff, block_per,link1, linkp, elev,
                H1, V1, H2, V2, dt, g, 0, np.sign(links1), [-1],
                friction, dVdx1, dVdx2, dVdt1, dVdt2)

    elif utype[0] == 'Pump':
        pumpc = pump[0]
        H[0], V[0] = pump_node(pumpc, link1, linkp,
                H1, V1, H2, V2, dt, g, 0,  np.sign(links1), [-1],
                friction, dVdx1, dVdx2, dVdt1, dVdt2)
    elif utype[0] == 'Valve':
        valvec = valve[0]
        H[0], V[0] = valve_node(valvec, link1, linkp,
                H1, V1, H2, V2, dt, g, 0,  np.sign(links1), [-1],
                friction, dVdx1, dVdx2, dVdt1, dVdt2)

    # Pipe end (outer boundary conditions )
    V1 = V0[n-1]; H1 = H0[n-1]
    dVdx1 = dVdx[n-1]
    dVdt1 = dVdt[n-1]
    if dtype[0] == 'Reservoir' or  dtype[0] == 'Tank':
        H[n], V[n] = rev_end (H1, V1, H[n], n, a, g, f, D, dt,
                            KD, friction, dVdx1, dVdt1)
    if  dtype[0] == 'Valve':
        H[n], V[n] = valve_end (H1, V1, V[n], n, a, g, f, D, dt,
                        KD, friction, dVdx1, dVdt1)
    if dtype[0] == 'Junction':
        elev = linkp.end_node.elevation
        H[n], V[n] = dead_end (linkp ,H1, V1, elev, n, a, g, f, D, dt,
                    KD, friction, dVdx1, dVdt1)


    return H, V