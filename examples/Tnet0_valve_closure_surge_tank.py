import tsnet
# Open an example network and create a transient model
inp_file = '/Users/luxing/Code/TSNet/examples/networks/Tnet0.inp'
tm = tsnet.network.TransientModel(inp_file)

# Set wavespeed
tm.set_wavespeed(1200.) # m/s
# Set time options
dt = 0.1
tf = 60   # simulation period [s]
tm.set_time(tf,dt)

# Set valve closure
tc = 0 # valve closure period [s]
ts = 2 # valve closure start time [s]
se = 0 # end open percentage [s]
m = 1 # closure constant [dimensionless]
valve_op = [tc,ts,se,m]
tm.valve_closure('3',valve_op)

# Initialize steady state simulation
t0 = 0. # initialize the simulation at 0 [s]
engine = 'DD' # demand driven simulator
tm = tsnet.simulation.Initializer(tm, t0, engine)

tank_height = 100  # tank height [m]
water_height = 50  # initial water level [m]
tank_node = '2'
tank_area = 100  # tank cross sectional area [m^2]
tm.add_surge_tank(tank_node, [tank_area], 'open')

# Transient simulation
results_obj = 'Tnet0' # name of the object for saving simulation results
tm1 = tsnet.simulation.MOCSimulator(tm, results_obj)


#%%
tm = tsnet.network.TransientModel(inp_file)

# Set wavespeed
tm.set_wavespeed(1200.) # m/s
# Set time options
tm.set_time(tf,dt)

# Set valve closure
tm.valve_closure('3',valve_op)

# Initialize steady state simulation
t0 = 0. # initialize the simulation at 0 [s]
engine = 'DD' # demand driven simulator
tm = tsnet.simulation.Initializer(tm, t0, engine)

tank_height = 100  # tank height [m]
water_height =40  # initial water level [m]
tank_node = '2'
tank_area = 10  # tank cross sectional area [m^2]
tm.add_surge_tank(tank_node, [tank_area,tank_height,water_height], 'closed')

# Transient simulation
results_obj = 'Tnet0' # name of the object for saving simulation results
tm2 = tsnet.simulation.MOCSimulator(tm, results_obj)

#%%
tm = tsnet.network.TransientModel(inp_file)

# Set wavespeed
tm.set_wavespeed(1200.) # m/s
# Set time options
tm.set_time(tf,dt)

# Set valve closure
tm.valve_closure('3',valve_op)

# Initialize steady state simulation
t0 = 0. # initialize the simulation at 0 [s]
engine = 'DD' # demand driven simulator
tm = tsnet.simulation.Initializer(tm, t0, engine)

# Transient simulation
results_obj = 'Tnet0' # name of the object for saving simulation results
tm3 = tsnet.simulation.MOCSimulator(tm, results_obj)


#%%
# report results
import matplotlib.pyplot as plt
node = '3'
norm_head_nl1 = tm1.get_node(node).head - tm1.get_node(node).head[0]
norm_head_nl2 = tm2.get_node(node).head - tm2.get_node(node).head[0]
norm_head_nl3 = tm3.get_node(node).head - tm3.get_node(node).head[0]

fig1 = plt.figure(figsize=(8,5), dpi=80, facecolor='w', edgecolor='k')
plt.plot(tm.simulation_timestamps,norm_head_nl1,'r-',label='w surge tank $A_s=10m^2$', linewidth=2.5)
plt.plot(tm.simulation_timestamps,norm_head_nl2,'g-',label='w surge tank $A_s=100m^2$', linewidth=2.5)
plt.plot(tm.simulation_timestamps,norm_head_nl3,'k-',label='wo surge tank', linewidth=2.5)
plt.xlim([tm.simulation_timestamps[0],tm.simulation_timestamps[-1]])
plt.title('Node %s' %node)
plt.xlabel("Time [s]")
plt.ylabel("Pressure Head Change [m]")
plt.legend(loc='best')
plt.grid(False)
plt.show()
# fig1.savefig('./docs/figures/tnet1_node.png', format='png',dpi=100)

#%%
# report results
import matplotlib.pyplot as plt
node1 = tm1.get_node(tank_node)
t1 = tm1.simulation_timestamps
node2 = tm2.get_node(tank_node)
t2 = tm1.simulation_timestamps
fig = plt.figure(figsize=(8,5), dpi=80, facecolor='w', edgecolor='k')
plt.plot(t1, node1.tank_flow_timeseries, 'r',label='w surge tank $A_s=10m^2$', linewidth=2.5)
plt.plot(t2, node2.tank_flow_timeseries, 'g',label='w surge tank $A_s=100m^2$', linewidth=2.5)
plt.xlim([t1[0],t1[-1]])
plt.xlabel("Time [s]")
plt.ylabel("Flow into tank [m^3s]")
plt.legend(loc='best')
plt.grid(True)
plt.show()
