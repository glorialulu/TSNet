#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 09:34:17 2020

@author: luxing
"""

# Define test-case options
initial_status_active = False
provide_curve = True
kl_fact = 0.1
#kl_fact = 0.001
minor_loss = 0
setting = 1
#minor_loss = 1e10
#setting = 1

# Initialize
import matplotlib
from matplotlib import pyplot as plt
import numpy
import pandas
import wntr
import tsnet
from tsnet.utils.valve_curve import valve_curve


def gpm_to_m3s(gpm):
    return gpm/(60*264.172)

def m3s_to_gpm(m3s):
    return m3s*(60*264.172)

def m_to_psi(m):
    return m*1.4219702063247

def psi_to_m(psi):
    return psi/1.4219702063247

# New network
wn = wntr.network.WaterNetworkModel()

# Add water source --> reservoir at atmos, with a pump to get significant pressure
wn.add_reservoir('res', base_head=0)
wn.add_junction('in')
wn.add_pipe('source', start_node_name='res', end_node_name='in', length=40, diameter=0.03175)
wn.add_junction('pumpjunc')
pumpcurve = []
pumpcurve.append((gpm_to_m3s(0),psi_to_m(60)))
pumpcurve.append((gpm_to_m3s(5),psi_to_m(45)))
pumpcurve.append((gpm_to_m3s(10),psi_to_m(25)))
wn.add_curve('pumpcurve', 'HEAD', pumpcurve)
wn.add_pump(name='supply', start_node_name='in', end_node_name='pumpjunc', pump_type='HEAD', pump_parameter='pumpcurve')

# Add reference point junction 'ref'
wn.add_junction('ref')
wn.add_pipe('to_ref', start_node_name='pumpjunc', end_node_name='ref', length=10, diameter=0.03175)

# Add valve
wn.add_junction('before_valve1')
wn.add_pipe('to_valve1', start_node_name='ref', end_node_name='before_valve1', length=25, diameter=0.03175)
wn.add_junction('after_valve1')
wn.add_valve('valve1', start_node_name='before_valve1', end_node_name='after_valve1', valve_type='FCV', minor_loss=minor_loss, setting=setting)
if initial_status_active:
    wn.get_link('valve1').initial_status = wntr.network.elements.LinkStatus.Active
else:
    wn.get_link('valve1').initial_status = wntr.network.elements.LinkStatus.Closed
wn.add_reservoir('atmos1', base_head=psi_to_m(0))
wn.add_pipe('out1', start_node_name='after_valve1', end_node_name='atmos1', length=15, diameter=0.03175/3)

# Prep plotting
fig = plt.figure(figsize=(8,5), dpi=80, facecolor='w', edgecolor='k')

# Build EPANET file
wn.write_inpfile('simp.inp')

# Construct transient model 
tm = tsnet.network.TransientModel('simp.inp')
tm.set_wavespeed(1200.) # m/s
tm.set_time(3) # simulation period [s]

# Extract default valve curve
po = numpy.linspace(100,0,11)
kl = valve_curve(po)

# Scale loss by requested factor and assemble the curve in the required format
kl = kl*kl_fact
curve = list(zip(po, kl))

# Confirm 'curve' contains the default curve
if kl_fact == 1:
    curve_test1 = valve_curve(po, coeff=(po, kl))
    curve_test2 = valve_curve(po, coeff=list(zip(*curve)))
    assert((curve_test1 == kl).all())
    assert((curve_test2 == kl).all())
    print('default valve_curve passed consistency check')

# Plot the valve curve
plt.clf()
plt.plot(po, kl, 'o')
plt.xlabel('% open')
plt.ylabel('loss')
plt.draw()
plt.savefig('valve_curve.png')

# Set valve opening 
tc = 0.2 # valve opening period [s]
ts = 1 # valve opening start time [s]
se = 1 # end open percentage [s]
m = 1 # open constant [dimensionless]
if provide_curve:
    tm.valve_opening('valve1',[tc,ts,se,m], curve=curve)
else:
    tm.valve_opening('valve1',[tc,ts,se,m])

# Initialize and run simulation
tm = tsnet.simulation.Initializer(tm, 0, 'DD')
tm = tsnet.simulation.MOCSimulator(tm, 'simpres', 'unsteady')

# Gather results
t = tm.simulation_timestamps
ref_flow_rate = tm.get_link('to_ref').end_node_flowrate

# Prep plotting
fig = plt.figure(figsize=(8,5), dpi=80, facecolor='w', edgecolor='k')

# Plot flow-rates
plt.clf()
plt.plot(t, m3s_to_gpm(ref_flow_rate), label='ref', linewidth=1.5)
plt.plot(t, m3s_to_gpm(tm.get_link('out1').end_node_flowrate), label='output1', linewidth=1.5)
plt.ylabel("flow-rate [gpm]")
plt.xlim([t[0],t[-1]])
plt.xlabel("Time [s]")
plt.legend(loc='best')
plt.draw()
plt.savefig('flowrate.png')

# Plot pressure
plt.clf()
plt.plot(t, m_to_psi(tm.get_node('ref').head), 'k', label='ref', linewidth=1.5)
plt.plot(t, m_to_psi(tm.get_node('after_valve1').head), label='output1', linewidth=1.5)
plt.ylabel("Pressure [PSI]")
plt.xlim([t[0],t[-1]])
plt.xlabel("Time [s]")
plt.legend(loc='best')
plt.draw()
plt.savefig('pressure.png')

print('Initial ref flow rate={}'.format(m3s_to_gpm(ref_flow_rate[0])))
print('Final ref flow rate={}'.format(m3s_to_gpm(ref_flow_rate[-1])))