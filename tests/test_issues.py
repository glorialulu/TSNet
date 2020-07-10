import matplotlib
# matplotlib.use('agg')
from matplotlib import pyplot as plt
import wntr
import tsnet


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

# Add components
wn.add_reservoir('res', base_head=0, coordinates=(-100,0))
wn.add_junction('in', coordinates=(-20,0))
wn.add_pipe('source', start_node_name='res', end_node_name='in', length=80, diameter=0.03175*2)
wn.add_junction('ref', coordinates=(0,0))
pumpcurve = []
pumpcurve.append((gpm_to_m3s(0),psi_to_m(39)))
pumpcurve.append((gpm_to_m3s(1),psi_to_m(30)))
pumpcurve.append((gpm_to_m3s(2),psi_to_m(20)))
wn.add_curve('pumpcurve', 'HEAD', pumpcurve)
wn.add_pump(name='supply', start_node_name='in', end_node_name='ref', pump_type='HEAD', pump_parameter='pumpcurve')
wn.add_reservoir('atmos', base_head=psi_to_m(0), coordinates=(40,0))
wn.add_junction('before_valve', coordinates=(20,0))
wn.add_pipe('to_valve', start_node_name='ref', end_node_name='before_valve', length=100, diameter=0.03175)
wn.add_junction('after_valve', base_demand=0, coordinates=(20,0))
wn.add_valve('valve', start_node_name='before_valve', end_node_name='after_valve', valve_type='FCV', setting=0)
wn.get_link('valve').status = wntr.network.elements.LinkStatus.Closed
wn.add_reservoir('atmos', base_head=0, coordinates=(40,0))
wn.add_pipe('out', start_node_name='after_valve', end_node_name='atmos', length=10, diameter=0.03175/3)

# Build EPANET file
wn.write_inpfile('simp.inp')

# Construct transient model 
tm = tsnet.network.TransientModel('simp.inp')
wn.get_link('valve').initial_status = wntr.network.elements.LinkStatus.Closed

# Set wavespeed
tm.set_wavespeed(1200.) # m/s

# Set time options
tf = 5   # simulation period [s]
tm.set_time(tf)

# Set valve opening 
tc = 0.1 # valve opening period [s]
ts = 0.5 # valve opening start time [s]
se = 0.1 # end open percentage [s]
m = 1 # open constant [dimensionless]
tm.valve_opening('valve',[tc,ts,se,m])

# Initialize steady state simulation
t0=0
tm = tsnet.simulation.Initializer(tm,t0,'DD')

# Transient simulation
results_obj = 'simpres'
friction = 'unsteady'
tm = tsnet.simulation.MOCSimulator(tm, results_obj, friction)

# Gather results
t = tm.simulation_timestamps
ref_head = tm.get_node('ref').head
output_head = tm.get_node('after_valve').head
ref_flow_rate = tm.get_link('to_valve').start_node_flowrate
output_flow_rate = tm.get_link('out').end_node_flowrate

# Prep plotting
fig = plt.figure(figsize=(8,5), dpi=80, facecolor='w', edgecolor='k')

# Plot flow-rates
plt.plot(t, m3s_to_gpm(ref_flow_rate), 'k', label='ref', linewidth=1.5)
plt.plot(t, m3s_to_gpm(output_flow_rate), 'r--', label='output', linewidth=1.5)
plt.ylabel("flow-rate [gpm]")
plt.xlim([t[0],t[-1]])
plt.xlabel("Time [s]")
plt.legend(loc='best')
plt.show()
plt.savefig('flowrate.png')
    
# Plot heads
plt.clf()
plt.plot(t, m_to_psi(ref_head), 'k', label='ref', linewidth=1.5)
plt.plot(t, m_to_psi(output_head), 'r--', label='output', linewidth=1.5)
plt.ylabel("Pressure [PSI]")
plt.xlim([t[0],t[-1]])
plt.xlabel("Time [s]")
plt.legend(loc='best')
plt.show()
plt.savefig('pressure.png')


wntr.graphics.plot_network(wn)
plt.show()



