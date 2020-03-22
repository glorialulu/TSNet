"""
The tsnet.simulation.main module contains function to perform
the workflow of read, discretize, initial, and transient
simulation for the given .inp file.

"""
from __future__ import print_function
from tsnet.network import  topology
from tsnet.simulation.single import inner_pipe, left_boundary, right_boundary
from tsnet.utils import valve_curve, memo, print_time_delta
from tsnet.utils import calc_parabola_vertex
import numpy as np
import warnings
from datetime import datetime
import pickle

@memo
def MOCSimulator(tm, results_obj='results', friction='steady'):
    """ MOC Main Function

    Parameters
    ----------
    tm : tsnet.network.model.TransientModel
        Network
    results_obj: string, optional
        the name of the results file, by default 'results'
    friction: string, optional
        friction model, e.g., 'steady', 'quasi-steady', 'unsteady',
        by default 'steady'
    Returns
    ------
    tm : tsnet.network.model.TransientModel
            Simulated network
    """

    # determine network topology
    links1, links2, utype, dtype = topology(tm)

    tt = ['x']
    tt.append(0)
    dt = tm.time_step
    tn = int(tm.simulation_period/tm.time_step)  # Total time steps
    # check whether input is legal
    if not friction in ['steady', 'unsteady', 'quasi-steady']:
        print ("Please specify a friction model from 'steady', 'unsteady', and 'quasi-steady'")

    # determine which node of the adjacent pipe should be call:
    # if the adjacent pipe is entering the junction, then -2
    # if the adjacent pipe is leaving the junction, then 1
    a = {1:-2, -1:1}
    b = {1:-1, -1:0}
    # generat a list of pipe
    p = []
    # results from last time step
    H = [0] * tm.num_pipes
    V = [0] * tm.num_pipes
    # results at current time step
    HN = [0] * tm.num_pipes
    VN = [0] * tm.num_pipes
    # results for local and convective
    #  instantaneous acceleration
    dVdt = [0] * tm.num_pipes
    dVdx = [0] * tm.num_pipes
    Hb = 10.3 # barometric head
    for _, pipe in tm.pipes():
        p.append(pipe)

    # initial condition
    for _, pipe in tm.pipes():
        pn = pipe.id-1
        H[pn] = pipe.initial_head
        V[pn] = pipe.initial_velocity
        if friction == 'unsteady':
            dVdt[pn] = np.zeros_like(V[pn])
            dVdx[pn] = np.diff(V[pn])/(pipe.length/pipe.number_of_segments)
        else:
            dVdt[pn] = np.zeros_like(V[pn])
            dVdx[pn] = np.zeros_like(V[pn][:-1])
    for _,node in tm.nodes():
        if node.pulse_status == True:
            node.base_demand_coeff = node.demand_coeff
        if node.transient_node_type == 'SurgeTank' or node.transient_node_type == 'Chamber':
            if node.transient_node_type == 'Chamber':
                m = 1.2
                Ha = node.initial_head - node.water_level + Hb # air pressure head
                Va = node.tank_shape[0]*(node.tank_height-node.water_level) # air volume
                node.air_constant = Ha * Va**m
                node.tank_shape.insert(2,node.air_constant)
            elif node.transient_node_type == 'SurgeTank':
                node.water_level = node.initial_head
                node.tank_shape.insert(1,node.water_level)
            node.water_level_timeseries = np.zeros(tn)
            node.tank_flow_timeseries = np.zeros(tn)
            node.water_level_timeseries[0] = node.water_level
    starttime = datetime.now()
    # Start Calculation
    for ts in range(1,tn):
        # check the discrepency between initial condition and the
        # first step in the transient simulation.
        if ts == 2:
            for _,pipe in tm.pipes():
                diff1 = pipe.start_node_head[1] - pipe.start_node_head[0]
                diff2 = pipe.end_node_head[1] - pipe.end_node_head[0]
                if abs(diff1)> 1e-1:
                    print('Initial condition discrepancy of pressure (%.4f m) on the %s node' %(diff1,pipe.start_node.name))
                if abs(diff2)> 1e-1:
                    print('Initial condition discrepancy of pressure (%.4f m) on the %s node'%(diff2,pipe.end_node.name))
        if ts == 3:
            timeperstep = (datetime.now() - starttime) /2.
            est = timeperstep *tn
            print ('Estimated simulation time %s' %est)

        t = ts*dt
        tt.append(t)
        tp = ts/tn*100
        if ts % int(tn/10) == 0 :
            print('Transient simulation completed %i %%...' %tp )
        # for burst node: emitter_coeff = burst_coeff[ts]
        for _,node in tm.nodes():
            if node.burst_status == True:
                node.emitter_coeff = node.burst_coeff[ts]
            if node.pulse_status == True:
                node.demand_coeff = node.base_demand_coeff*(1.+node.pulse_coeff[ts])

        # initialize the results at this time step
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
                    pump[0] = [tm.links[utype[pn][1]].curve_coef, "d"]
                    if pipe.start_node.name == tm.links[utype[pn][1]].start_node.name:
                        pump[0][1] = "s" # suction side
                    # calculate the coordinate of the three points
                    # based on the pump speed
                    if tm.links[utype[pn][1]].operating == True:
                        points = tm.links[utype[pn][1]].get_pump_curve().points
                        po = tm.links[utype[pn][1]].operation_rule[ts]
                        points=[(i*po,j*po**2) for (i,j) in points]
                        pump[0][0] = calc_parabola_vertex(points)

                elif utype[pn][0] == 'Valve':
                    # determine valve friction coefficients based on
                    # open percentage
                    if tm.links[utype[pn][1]].operating == True:
                        valve[0] = valve_curve(tm.links[utype[pn][1]].operation_rule[ts]*100,
                        tm.links[utype[pn][1]].valve_coeff)
                    else :
                         valve[0] = valve_curve(100,tm.links[utype[pn][1]].valve_coeff)
                # downstream
                if dtype[pn][0] == 'Pump':
                    pump[1] = [tm.links[dtype[pn][1]].curve_coef,"d"]
                    if pipe.end_node.name == tm.links[dtype[pn][1]].start_node.name:
                        pump[1][1] = "s" # suction side
                    if tm.links[dtype[pn][1]].operating == True:
                        points = tm.links[dtype[pn][1]].get_pump_curve().points
                        po = tm.links[dtype[pn][1]].operation_rule[ts]
                        points=[(i*po,j*po**2) for (i,j) in points]
                        pump[1][0] = calc_parabola_vertex(points)

                elif dtype[pn][0] == 'Valve':
                    if tm.links[dtype[pn][1]].operating == True:
                        valve[1] = valve_curve(tm.links[dtype[pn][1]].operation_rule[ts]*100,
                        tm.links[dtype[pn][1]].valve_coeff)
                    else :
                         valve[1] = valve_curve(100,tm.links[dtype[pn][1]].valve_coeff)

                HN[pn], VN[pn] = inner_pipe(pipe, pn, dt,
                     links1[pn], links2[pn], utype[pn], dtype[pn], p,
                     H[pn], V[pn], HN[pn], VN[pn],
                     [H[abs(i)-1][a[np.sign(i)]] for i in links1[pn]],
                     [V[abs(i)-1][a[np.sign(i)]] for i in links1[pn]],
                     [H[abs(i)-1][a[np.sign(i)]] for i in links2[pn]],
                     [V[abs(i)-1][a[np.sign(i)]] for i in links2[pn]],
                     pump, valve, friction, dVdt[pn], dVdx[pn],
                     [dVdt[abs(i)-1][a[np.sign(i)]] for i in links1[pn]],
                     [dVdx[abs(i)-1][b[np.sign(i)]] for i in links1[pn]],
                     [dVdt[abs(i)-1][a[np.sign(i)]] for i in links2[pn]],
                     [dVdx[abs(i)-1][b[np.sign(i)]] for i in links2[pn]])
                # record results
                pipe.start_node_velocity[ts] = VN[pn][0]
                pipe.end_node_velocity[ts] = VN[pn][-1]
                pipe.start_node_flowrate[ts] = VN[pn][0]*pipe.area
                pipe.end_node_flowrate[ts] = VN[pn][-1]*pipe.area
                pipe.start_node_head[ts] = HN[pn][0]
                pipe.end_node_head[ts] = HN[pn][-1]

                if pipe.start_node.transient_node_type == 'Junction':
                    if HN[pn][0] - pipe.start_node.elevation >0:
                        h = HN[pn][0] - pipe.start_node.elevation
                        pipe.start_node.demand_discharge[ts] = pipe.start_node.demand_coeff * np.sqrt(h)
                        pipe.start_node.emitter_discharge[ts] = pipe.start_node.emitter_coeff * np.sqrt(h)
                    else: # assume reverse flow preventer installed
                        pipe.start_node.emitter_discharge[ts] = 0.
                        pipe.start_node.demand_discharge[ts] = 0.
                        warnings.warn("Negative pressure on node %s.\
                        Backflow stopped by reverse flow preventer." %pipe.start_node.name)

                if pipe.end_node.transient_node_type == 'Junction':
                    if HN[pn][-1]  -pipe.end_node.elevation >0:
                        h = HN[pn][-1] -pipe.end_node.elevation
                        pipe.end_node.emitter_discharge[ts] = pipe.end_node.emitter_coeff * np.sqrt(h)
                        pipe.end_node.demand_discharge[ts] = pipe.end_node.demand_coeff * np.sqrt(h)
                    else: # assume reverse flow preventer installed
                        pipe.end_node.emitter_discharge[ts] = 0.
                        pipe.end_node.demand_discharge[ts] = 0.
                        warnings.warn("Negative pressure on node %s.\
                            Backflow stopped by reverse flow preventer." %pipe.start_node.name)

            # left boundary pipe
            elif not links1[pn] or links1[pn] == ['End']:
                pump = [[],[]]; valve = [0,0]
                # LEFT BOUNDARY
                if utype[pn][0] == 'Reservoir' or utype[pn][0] == 'Tank':
                    # head B.C.
                    HN[pn][0] = pipe.initial_head[0]
                elif utype[pn][0] == 'Junction':
                    VN[pn][0] = pipe.initial_velocity[0]
                elif utype[pn][0] == 'Valve':
                    if tm.links[utype[pn][1]].operating == True:
                        # velocity B.C.
                        VN[pn][0] = pipe.initial_velocity[0] * \
                            tm.links[utype[pn][1]].operation_rule[ts]
                    else :
                        VN[pn][0] = pipe.initial_velocity[0]
                elif utype[pn][0] == 'Pump':
                    # source pump
                    # pump[0][0]: elevation of the reservoir/tank
                    # pump[0][1]: three points for pump characteristic curve
                    pump[0] = [[tm.links[utype[pn][1]].start_node.initial_head][0],
                         tm.links[utype[pn][1]].curve_coef]
                    if tm.links[utype[pn][1]].operating == True:
                        points = tm.links[utype[pn][1]].get_pump_curve().points
                        po = tm.links[utype[pn][1]].operation_rule[ts]
                        points= [(i*po,j*po**2) for (i,j) in points]
                        pump[0][1] = calc_parabola_vertex(points)
                else:
                     warnings.warn ('Pipe %s miss %s upstream.' %(pipe, utype[pn][0]))

                # RIGHT BOUNDARY
                if dtype[pn][0] == 'Pump':
                    pump[1] = [tm.links[dtype[pn][1]].curve_coef,"d"]
                    if pipe.end_node.name == tm.links[dtype[pn][1]].start_node.name:
                        pump[1][1] = "s" # suction side
                    if tm.links[dtype[pn][1]].operating == True:
                        points = tm.links[dtype[pn][1]].get_pump_curve().points
                        po = tm.links[dtype[pn][1]].operation_rule[ts]
                        points=[(i*po,j*po**2) for (i,j) in points]
                        pump[1][0] = calc_parabola_vertex(points)

                elif dtype[pn][0] == 'Valve':
                    if tm.links[dtype[pn][1]].operating == True:
                        valve[1] = valve_curve(tm.links[dtype[pn][1]].operation_rule*100,
                        tm.links[dtype[pn][1]].valve_coeff)
                    else :
                         valve[1] = valve_curve(100, tm.links[dtype[pn][1]].valve_coeff)
                    # if also the right valve end
                    if links2[pn] == ['End']:
                        links2[pn] = []

                elif dtype[pn][0] == 'Junction':
                    VN[pn][-1] = pipe.initial_velocity[-1]

                HN[pn], VN[pn] = left_boundary(pipe, pn,
                      HN[pn], VN[pn], H[pn], V[pn],
                     links2[pn], p, pump, valve, dt,
                     [H[abs(i)-1][a[np.sign(i)]] for i in links2[pn]],
                     [V[abs(i)-1][a[np.sign(i)]] for i in links2[pn]],
                     utype[pn], dtype[pn],
                     friction, dVdt[pn], dVdx[pn],
                     [dVdt[abs(i)-1][a[np.sign(i)]] for i in links2[pn]],
                     [dVdx[abs(i)-1][b[np.sign(i)]] for i in links2[pn]],)
                # record results
                pipe.start_node_velocity[ts] = VN[pn][0]
                pipe.end_node_velocity[ts] = VN[pn][-1]
                pipe.start_node_head[ts] = HN[pn][0]
                pipe.end_node_head[ts] = HN[pn][-1]
                pipe.start_node_flowrate[ts] = VN[pn][0]*pipe.area
                pipe.end_node_flowrate[ts] = VN[pn][-1]*pipe.area

                try:
                    if HN[pn][0]- pipe.start_node.elevation >0:
                        h = HN[pn][0]- pipe.start_node.elevation
                        pipe.start_node.demand_discharge[ts] = pipe.start_node.demand_coeff * np.sqrt(h)
                        pipe.start_node.emitter_discharge[ts] = pipe.start_node.emitter_coeff * np.sqrt(h)
                    else: # assume reverse flow preventer installed
                        pipe.start_node.emitter_discharge[ts] = 0.
                        pipe.start_node.demand_discharge[ts] = 0.
                        warnings.warn("Negative pressure on node %s.\
                        Backflow stopped by reverse flow preventer." %pipe.start_node.name)
                except:
                    pass

                try:
                    if HN[pn][-1]-pipe.end_node.elevation >0:
                        h = HN[pn][-1]-pipe.end_node.elevation
                        pipe.end_node.emitter_discharge[ts] = pipe.end_node.emitter_coeff * np.sqrt(h)
                        pipe.end_node.demand_discharge[ts] = pipe.end_node.demand_coeff * np.sqrt(h)
                    else: # assume reverse flow preventer installed
                        pipe.end_node.emitter_discharge[ts] = 0.
                        pipe.end_node.demand_discharge[ts] = 0.
                        warnings.warn("Negative pressure on node %s.\
                            Backflow stopped by reverse flow preventer." %pipe.start_node.name)
                except:
                    pass

            #  right boundary pipe
            elif not links2[pn] or links2[pn] == ['End']:
                pump = [[],[]]; valve = [0,0]
                # RIGHT boundary
                if dtype[pn][0] == 'Reservoir' or dtype[pn][0] == 'Tank':
                    HN[pn][-1]   =  pipe.initial_head[-1] # head of reservoir
                elif dtype[pn][0] == 'Junction':
                    VN[pn][-1] = pipe.initial_velocity[-1]
                elif dtype[pn][0] == 'Valve':
                    if tm.links[dtype[pn][1]].operating == True:
                        # valve velocity condition
                        VN[pn][-1] = pipe.initial_velocity[-1]* \
                        tm.links[dtype[pn][1]].operation_rule[ts]
                    else :
                        VN[pn][-1] = pipe.initial_velocity[-1]
                # source pump
                elif dtype[pn][0] == 'Pump':
                    # pump[1][0]: elevation of the reservoir/tank
                    # pump[1][1]: three points for pump characteristic curve
                    pump[1] = [[tm.links[utype[pn][1]].end_node.initial_head][0],
                         tm.links[dtype[pn][1]].curve_coef]
                    if tm.links[dtype[pn][1]].operating == True:
                        points = tm.links[dtype[pn][1]].get_pump_curve().points
                        po = tm.links[dtype[pn][1]].operation_rule[ts]
                        points=[(i*po,j*po**2) for (i,j) in points]
                        pump[1][1] = calc_parabola_vertex(points)
                else :
                     warnings.warn('Pipe %s miss %s downstream.' %(pipe, dtype[pn][0]))
                # LEFT boundary
                if utype[pn][0] == 'Pump':
                    pump[0] = [tm.links[utype[pn][1]].curve_coef,"d"]
                    if pipe.start_node.name == tm.links[utype[pn][1]].start_node.name:
                        pump[0][1] = "s" # suction side
                    if tm.links[utype[pn][1]].operating == True:
                        points = tm.links[utype[pn][1]].get_pump_curve().points
                        po = tm.links[utype[pn][1]].operation_rule[ts]
                        points=[(i*po,j*po**2) for (i,j) in points]
                        pump[0][0] = calc_parabola_vertex(points)

                elif utype[pn][0] == 'Valve':
                    if tm.links[utype[pn][1]].operating == True:
                        valve[0] = valve_curve(tm.links[utype[pn][1]].operation_rule[ts]*100,
                        tm.links[utype[pn][1]].valve_coef)
                    else :
                         valve[0] = tm.links[utype[pn][1]].valve_curve(100,
                         tm.links[utype[pn][1]].valve_coef)

                HN[pn], VN[pn] = right_boundary(pipe, pn,
                     H[pn], V[pn], HN[pn], VN[pn],
                     links1[pn], p, pump, valve,  dt,
                     [H[abs(i)-1][a[np.sign(i)]] for i in links1[pn]],
                     [V[abs(i)-1][a[np.sign(i)]] for i in links1[pn]],
                     utype[pn], dtype[pn],
                     friction, dVdt[pn], dVdx[pn],
                     [dVdt[abs(i)-1][a[np.sign(i)]] for i in links1[pn]],
                     [dVdx[abs(i)-1][b[np.sign(i)]] for i in links1[pn]],)
                # record results
                pipe.start_node_velocity[ts] = VN[pn][0]
                pipe.end_node_velocity[ts] = VN[pn][-1]
                pipe.start_node_head[ts] = HN[pn][0]
                pipe.end_node_head[ts] = HN[pn][-1]
                pipe.start_node_flowrate[ts] = VN[pn][0]*pipe.area
                pipe.end_node_flowrate[ts] = VN[pn][-1]*pipe.area

                try:
                    if HN[pn][0]- pipe.start_node.elevation >0:
                        h = HN[pn][0]- pipe.start_node.elevation
                        pipe.start_node.demand_discharge[ts] = pipe.start_node.demand_coeff * np.sqrt(h)
                        pipe.start_node.emitter_discharge[ts] = pipe.start_node.emitter_coeff * np.sqrt(h)
                    else: # assume reverse flow preventer installed
                        pipe.start_node.emitter_discharge[ts] = 0.
                        pipe.start_node.demand_discharge[ts] = 0.
                        warnings.warn("Negative pressure on node %s.\
                        Backflow stopped by reverse flow preventer." %pipe.start_node.name)
                except:
                    pass

                try:
                    if HN[pn][-1]-pipe.end_node.elevation >0:
                        h = HN[pn][-1]-pipe.end_node.elevation
                        pipe.end_node.emitter_discharge[ts] = pipe.end_node.emitter_coeff * np.sqrt(h)
                        pipe.end_node.demand_discharge[ts] = pipe.end_node.demand_coeff * np.sqrt(h)
                    else: # assume reverse flow preventer installed
                        pipe.end_node.emitter_discharge[ts] = 0.
                        pipe.end_node.demand_discharge[ts] = 0.
                        warnings.warn("Negative pressure on node %s.\
                            Backflow stopped by reverse flow preventer." %pipe.start_node.name)
                except:
                    pass

        # march in time
        for _, pipe in tm.pipes():
            pn = pipe.id-1
            # calculate instantaneous local acceleration
            # only for unsteady friction factor
            if friction == 'unsteady':
                dVdt[pn] = (VN[pn] - V[pn] )/dt
                dVdx[pn] =  np.diff(V[pn])/(pipe.length/pipe.number_of_segments)
            H[pn] = HN[pn]
            V[pn] = VN[pn]

        for _,node in tm.nodes():
            if node.transient_node_type == 'SurgeTank' or node.transient_node_type == 'Chamber':
                node.tank_shape[-2] = max(node.water_level,0)
                node.tank_shape[-1] = node.tank_flow
                node.water_level_timeseries[ts] = max(node.water_level,0)
                node.tank_flow_timeseries[ts] = node.tank_flow

    for _,pipe in tm.pipes():
        if not isinstance(pipe.start_node.head, np.ndarray):
            pipe.start_node.head = np.copy(pipe.start_node_head)
        if not isinstance(pipe.end_node.head, np.ndarray):
            pipe.end_node.head = np.copy(pipe.end_node_head)

    tm.simulation_timestamps = tt[1:]

    # save object to file
    import pickle
    filehandler = open(results_obj +'.obj','wb')
    pickle.dump(tm, filehandler)

    """TO read:
    import pickle
    file = open('results.obj', 'rb')
    tm = pickle.load(file)
    """
    return tm
