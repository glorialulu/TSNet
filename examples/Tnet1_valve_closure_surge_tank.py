import tsnet
import numpy as np 
# Open an example network and create a transient model
inp_file = 'networks/Tnet1.inp'
# Set valve closure
tc = 0.6 # valve closure period [s]
ts = 2  # valve closure start time [s]
se = 0 # end open percentage [s]
m = 1 # closure constant [dimensionless]
valve_op = [tc,ts,se,m]
#%%
tm = tsnet.network.TransientModel(inp_file)

# Set wavespeed
tm.set_wavespeed(1200.) # m/s
# Set time options
dt = 0.1  # time step [s], if not given, use the maximum allowed dt
tf = 60   # simulation period [s]
tm.set_time(tf,dt)

# Set valve closure
tm.valve_closure('VALVE', valve_op)

# Initialize steady state simulation
t0 = 0. # initialize the simulation at 0 [s]
engine = 'DD' # demand driven simulator
tm = tsnet.simulation.Initializer(tm, t0, engine)
tank_height = 10  # tank height [m]
water_height = 5  # initial water level [m]
tank_node = 'N4'
tank_area = 10   # tank cross sectional area [m^2]
tm.add_surge_tank(tank_node, [tank_area,tank_height,water_height], 'closed')
# Transient simulation
results_obj = 'Tnet1' # name of the object for saving simulation results
friction = 'steady'
tm1 = tsnet.simulation.MOCSimulator(tm, results_obj,friction)

#%%
tm = tsnet.network.TransientModel(inp_file)
# Set wavespeed
tm.set_wavespeed(1200.) # m/s
# Set time options
tm.set_time(tf,dt)

tm.valve_closure('VALVE', valve_op)

# Initialize steady state simulation
t0 = 0. # initialize the simulation at 0 [s]
engine = 'DD' # demand driven simulator
tm = tsnet.simulation.Initializer(tm, t0, engine)
# Transient simulation
results_obj = 'Tnet1' # name of the object for saving simulation results
friction = 'steady'
tm2= tsnet.simulation.MOCSimulator(tm, results_obj,friction)
#%%
# report results
import matplotlib.pyplot as plt
node = 'N2'
head1 = tm1.get_node(node).head
t1 = tm1.simulation_timestamps
head2 = tm2.get_node(node).head
t2 = tm2.simulation_timestamps
fig = plt.figure(figsize=(8,5), dpi=80, facecolor='w', edgecolor='k')
plt.plot(t1, head1, 'r',label='w surge tank', linewidth=2.5)
plt.plot(t2, head2, 'k', label='wo surge tank',  linewidth=2.5)
plt.xlim([t1[0],t1[-1]])
plt.title(node)
plt.xlabel("Time [s]")
plt.ylabel("Pressure Head [m] ")
plt.legend(loc='best')
plt.show()

# %%
# report results
import matplotlib.pyplot as plt
node = tm1.get_node(tank_node)
t1 = tm1.simulation_timestamps
fig = plt.figure(figsize=(8,5), dpi=80, facecolor='w', edgecolor='k')
plt.plot(t1, node.tank_flow_timeseries, 'k',label='w surge tank', linewidth=2.5)
plt.xlim([t1[0],t1[-1]])
plt.xlabel("Time [s]")
plt.ylabel("Flow into tank [m^3/s]")
plt.legend(loc='best')
plt.show()

fig = plt.figure(figsize=(8,5), dpi=80, facecolor='w', edgecolor='k')
plt.plot(t1, node.water_level_timeseries, 'k',label='w surge tank', linewidth=2.5)
plt.xlim([t1[0],t1[-1]])
plt.xlabel("Time [s]")
plt.ylabel("Water level in tank [m]")
plt.legend(loc='best')
plt.show()
