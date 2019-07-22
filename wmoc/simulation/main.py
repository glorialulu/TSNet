"""
The wmoc.simulation.main module contains function to perform
the workflow of read, discretize, initial, and transient
simulation for the given .inp file.

"""
from __future__ import print_function
from wmoc.network import  topology
from wmoc.simulation.single import inner_pipe, left_boundary, right_boundary
from wmoc.utils import valve_curve, memo
import numpy as np
import warnings

@memo
def MOCSimulator(tm):
    """ MOC Main Function

    Parameters
    ----------
    tm : wmoc.network.model.TransientModel
        Network

    Returns
    ------
    tm : wmoc.network.model.TransientModel
            Simulated network
    """
    # determine network topology
    links1, links2, utype, dtype = topology(tm)

    tt = ['x']
    tt.append(0)
    dt = tm.time_step
    tn = int(tm.simulation_peroid/tm.time_step)  # Total time steps

    # determine which node of the adjacant pipe should be call:
    # if the adjacant pipe is entering the junction, then -2
    # if the adjacent pipe is leaving the junction, then 1
    a = {1:-2, -1:1}
    # generat a list of pipe
    p = []
    # results from last time step
    H = [0] * tm.num_pipes
    V = [0] * tm.num_pipes
    # results at current time step
    HN = [0] * tm.num_pipes
    VN = [0] * tm.num_pipes

    for _, pipe in tm.pipes():
        p.append(pipe)

    # initial condition
    for _, pipe in tm.pipes():
        pn = pipe.id-1
        H[pn] = pipe.initial_head
        V[pn] = pipe.initial_velocity

    # Start Claculation
    for ts in range(1,tn):
        t = ts*dt
        tt.append(t)

        # for burst node: emitter_coeff = burst_coeff[ts]
        for _,node in tm.nodes():
            if node.burst_status == True:
                node.emitter_coeff = node.burst_coeff[ts]

        # initilaize the results at this time step
        for _, pipe in tm.pipes():
            pn = pipe.id-1
            HN[pn] =  np.zeros_like(H[pn])
            VN[pn] =  np.zeros_like(V[pn])

        for _,pipe in tm.pipes():
            pn = pipe.id-1
            # Assumption:
            # when a pipe is connected with a pump or valve,
            # the connection is not branch junction.

            # inner pipes
            if links1[pn] and links2[pn] and \
                links1[pn] != ['End'] and links2[pn] != ['End']:
                # list to store information about pump and vale
                # pump[0] and valve[0] for upstream elemnets
                # pump[1] and valve[1] for downstream elements
                pump = [[],[]]; valve = [0,0]
                # upstream
                if utype[pn][0] == 'Pump':
                    # three points for pump charatersitics curve
                    pump[0] = tm.links[utype[pn][1]].get_pump_curve().points
                    # calculate the coordinate of the three points
                    # based on the pump speed
                    if tm.links[utype[pn][1]].operating == True:
                        po = tm.links[utype[pn][1]].operation_rule[ts]
                        pump[0]=[(i*po,j*po**2) for (i,j) in pump[0]]

                elif utype[pn][0] == 'Valve':
                    # determine valve fricton coefficients based on
                    # open percentage
                    if tm.links[utype[pn][1]].operating == True:
                        valve[0] = valve_curve(tm.links[utype[pn][1]].operation_rule[ts]*100)
                    else :
                         valve[0] = valve_curve(100)
                # downstream
                if dtype[pn][0] == 'Pump':
                    pump[1] = tm.links[dtype[pn][1]].get_pump_curve().points
                    if tm.links[dtype[pn][1]].operating == True:
                        po = tm.links[dtype[pn][1]].operation_rule[ts]
                        pump[1]=[(i*po,j*po**2) for (i,j) in pump[1]]

                elif dtype[pn][0] == 'Valve':
                    if tm.links[dtype[pn][1]].operating == True:
                        valve[1] = valve_curve(tm.links[dtype[pn][1]].operation_rule[ts]*100)
                    else :
                         valve[1] = valve_curve(100)

                HN[pn], VN[pn] = inner_pipe(pipe, pn, dt,
                     links1[pn], links2[pn], utype[pn], dtype[pn], p,
                     H[pn], V[pn], HN[pn], VN[pn],
                     [H[abs(i)-1][a[np.sign(i)]] for i in links1[pn]],
                     [V[abs(i)-1][a[np.sign(i)]] for i in links1[pn]],
                     [H[abs(i)-1][a[np.sign(i)]] for i in links2[pn]],
                     [V[abs(i)-1][a[np.sign(i)]] for i in links2[pn]],
                     pump, valve)
                # record results
                pipe.start_node_velocity[ts] = VN[pn][0]
                pipe.end_node_velocity[ts] = VN[pn][-1]
                pipe.start_node_head[ts] = HN[pn][0]
                pipe.end_node_head[ts] = HN[pn][-1]

                # pipe.start_node.head[ts] = HN[pn][0]
                # pipe.end_node.head[ts] = HN[pn][-1]
                pipe.start_node.demand_discharge[ts] = pipe.start_node.demand_coeff * np.sqrt(HN[pn][0])
                pipe.end_node.demand_discharge[ts] =  pipe.end_node.demand_coeff * np.sqrt(HN[pn][-1])
                pipe.start_node.emitter_discharge[ts] = pipe.start_node.emitter_coeff * np.sqrt(HN[pn][0])
                pipe.end_node.emitter_discharge[ts] = pipe.end_node.emitter_coeff * np.sqrt(HN[pn][-1])

            # left boundary pipe
            elif not links1[pn] or links1[pn] == ['End']:
                pump = [[],[]]; valve = [0,0]
                # LEFT BOUNDARY
                if utype[pn][0] == 'Reservoir':
                    # head B.C.
                    HN[pn][0]   =  tm.nodes[utype[pn][1]].base_head
                elif utype[pn][0] == 'Tank':
                    # head B.C.
                    HN[pn][0]   =  tm.nodes[utype[pn][1]].head
                elif utype[pn][0] == 'Junction':
                    VN[pn][0] = pipe.initial_velocity[0]
                elif utype[pn][0] == 'Valve':
                    if tm.links[utype[pn][1]].operating == True:
                        # velocity B.C.
                        VN[pn][0]   = pipe.initial_velocity[0] * \
                            tm.links[utype[pn][1]].operation_rule[ts]
                    else :
                        VN[pn][0]   = pipe.initial_velocity[0]
                # source pump
                elif utype[pn][0] == 'Pump':
                    # pump[0][0]: elevation of the reservoir/tank
                    # pump[0][1]: three points for pump characteristic curve
                    pump[0] = [tm.links[utype[pn][1]].start_node.base_head,
                         tm.links[utype[pn][1]].get_pump_curve().points]

                    if tm.links[utype[pn][1]].operating == True:
                        po = tm.links[utype[pn][1]].operation_rule[ts]
                        pump[0][1]=[(i*po,j*po**2) for (i,j) in pump[0][1]]
                        print('Operating pump %s' %utype[pn][1])
                else:
                     warnings.warn ('Pipe %s miss %s upstream.' %(pipe, utype[pn][0]))

                # RIGHT BOUNDARY
                if dtype[pn][0] == 'Pump':
                    pump[1] = tm.links[dtype[pn][1]].get_pump_curve().points
                    if tm.links[dtype[pn][1]].operating == True:
                        po = tm.links[dtype[pn][1]].operation_rule[ts]
                        pump[1]=[(i*po[ts],j*po[ts]**2) for (i,j) in pump[1]]

                elif dtype[pn][0] == 'Valve':
                    if tm.links[dtype[pn][1]].operating == True:
                        valve[1] = valve_curve(tm.links[dtype[pn][1]].operation_rule*100)
                    else :
                         valve[1] = valve_curve(100)

                HN[pn], VN[pn] = left_boundary(pipe, pn,
                      HN[pn], VN[pn], H[pn], V[pn],
                     links2[pn], p, pump, valve, dt,
                     [H[abs(i)-1][a[np.sign(i)]] for i in links2[pn]],
                     [V[abs(i)-1][a[np.sign(i)]] for i in links2[pn]],
                     utype[pn], dtype[pn])
                # record results
                pipe.start_node_velocity[ts] = VN[pn][0]
                pipe.end_node_velocity[ts] = VN[pn][-1]
                pipe.start_node_head[ts] = HN[pn][0]
                pipe.end_node_head[ts] = HN[pn][-1]

                # pipe.start_node.head[ts] = HN[pn][0]
                # pipe.end_node.head[ts] = HN[pn][-1]
                pipe.start_node.demand_discharge[ts] = pipe.start_node.demand_coeff * np.sqrt(HN[pn][0])
                pipe.end_node.demand_discharge[ts] =  pipe.end_node.demand_coeff * np.sqrt(HN[pn][-1])
                pipe.start_node.emitter_discharge[ts] = pipe.start_node.emitter_coeff * np.sqrt(HN[pn][0])
                pipe.end_node.emitter_discharge[ts] = pipe.end_node.emitter_coeff * np.sqrt(HN[pn][-1])

            #  right boundary pipe
            elif not links2[pn] or links2[pn] == ['End']:
                pump = [[],[]]; valve = [0,0]
                # RIGHT boundary
                if dtype[pn][0] == 'Reservoir':
                    HN[pn][-1]   =  tm.nodes[dtype[pn][1]].base_head # head of reservoir
                elif dtype[pn][0] == 'Tank':
                    HN[pn][-1]   =  tm.nodes[dtype[pn][1]].head # head of tank
                elif dtype[pn][0] == 'Junction':
                    VN[pn][-1] = pipe.initial_velocity[-1]
                elif dtype[pn][0] == 'Valve':
                    if tm.links[dtype[pn][1]].operating == True:
                        # valve velocity condition
                        VN[pn][-1] = pipe.initial_velocity[-1]* \
                        tm.links[dtype[pn][1]].operation_rule[ts]
                    else :
                        VN[pn][-1] = pipe.initial_velocity[-1]
                else :
                     warnings.warn('Pipe %s miss %s downstream.' %(pipe, dtype[pn][0]))
                # LEFT boundary
                # source pump
                if utype[pn][0] == 'Pump':
                    # pump[0][0]: elevation of the reservoir/tank
                    # pump[0][1]: three points for pump characteristic curve
                    pump[0] = tm.links[utype[pn][1]].get_pump_curve().points
                    if tm.links[utype[pn][1]].operating == True:
                        po = tm.links[utype[pn][1]].operation_rule[ts]
                        pump[0]=[(i*po,j*po**2) for (i,j) in pump[0]]

                elif utype[pn][0] == 'Valve':
                    if tm.links[utype[pn][1]].operating == True:
                        valve[0] = valve_curve(tm.links[utype[pn][1]].operation_rule[ts]*100)
                    else :
                         valve[0] = valve_curve(100)

                HN[pn], VN[pn] = right_boundary(pipe, pn,
                     H[pn], V[pn], HN[pn], VN[pn],
                     links1[pn], p, pump, valve,  dt,
                     [H[abs(i)-1][a[np.sign(i)]] for i in links1[pn]],
                     [V[abs(i)-1][a[np.sign(i)]] for i in links1[pn]],
                     utype[pn], dtype[pn])
                # record results
                pipe.start_node_velocity[ts] = VN[pn][0]
                pipe.end_node_velocity[ts] = VN[pn][-1]
                pipe.start_node_head[ts] = HN[pn][0]
                pipe.end_node_head[ts] = HN[pn][-1]

                # pipe.start_node.head[ts] = HN[pn][0]
                # pipe.end_node.head[ts] = HN[pn][-1]
                pipe.start_node.demand_discharge[ts] = pipe.start_node.demand_coeff * np.sqrt(HN[pn][0])
                pipe.end_node.demand_discharge[ts] =  pipe.end_node.demand_coeff * np.sqrt(HN[pn][-1])
                pipe.start_node.emitter_discharge[ts] = pipe.start_node.emitter_coeff * np.sqrt(HN[pn][0])
                pipe.end_node.emitter_discharge[ts] = pipe.end_node.emitter_coeff * np.sqrt(HN[pn][-1])

        # march in time
        for _, pipe in tm.pipes():
            pn = pipe.id-1
            H[pn] = HN[pn]
            V[pn] = VN[pn]

    for _,pipe in tm.pipes():
        if not isinstance(pipe.start_node.head, np.ndarray ):
            pipe.start_node.head = np.copy(pipe.start_node_head)
        if not isinstance(pipe.end_node.head, np.ndarray):
            pipe.end_node.head = np.copy(pipe.end_node_head)

    tm.simulation_timestamps = tt[1:]
    return tm
