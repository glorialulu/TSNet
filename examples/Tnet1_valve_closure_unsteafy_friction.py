import tsnet
import numpy as np 
# Open an example network and create a transient model
inp_file = 'networks/Tnet1.inp'
tm = tsnet.network.TransientModel(inp_file)

# Set wavespeed
tm.set_wavespeed(1200.) # m/s
# Set time options
dt = 0.02  # time step [s], if not given, use the maximum allowed dt
tf = 60   # simulation period [s]
tm.set_time(tf,dt)

# Set valve closure
tc = 0.6 # valve closure period [s]
ts = 0 # valve closure start time [s]
se = 0 # end open percentage [s]
m = 1 # closure constant [dimensionless]
valve_op = [tc,ts,se,m]
percent_open = np.linspace(100,0,11)
kl = [1/0.2, 2.50, 1.25, 0.625, 0.333, 0.17,
            0.100, 0.0556, 0.0313, 0.0167, 0.0]
curve = [(percent_open[i], kl[i]) for i in range(len(kl))]
tm.valve_closure('VALVE', valve_op, curve)

# Initialize steady state simulation
t0 = 0. # initialize the simulation at 0 [s]
engine = 'DD' # demand driven simulator
tm = tsnet.simulation.Initializer(tm, t0, engine)

# Transient simulation
results_obj = 'Tnet1' # name of the object for saving simulation results
friction = 'steady'
tm1 = tsnet.simulation.MOCSimulator(tm, results_obj,friction)


#%%
tm = tsnet.network.TransientModel(inp_file)

# Set wavespeed
tm.set_wavespeed(1200.) # m/s
# Set time options
tf = 60   # simulation period [s]
tm.set_time(tf,dt)

# Set valve closure
tc = 0.6 # valve closure period [s]
ts = 0 # valve closure start time [s]
se = 0 # end open percentage [s]
m = 1 # closure constant [dimensionless]
valve_op = [tc,ts,se,m]
tm.valve_closure('VALVE',valve_op)

# Initialize steady state simulation
t0 = 0. # initialize the simulation at 0 [s]
engine = 'DD' # demand driven simulator
tm = tsnet.simulation.Initializer(tm, t0, engine)

# Transient simulation
results_obj = 'Tnet1' # name of the object for saving simulation results
friction = 'quasi-steady'
tm2 = tsnet.simulation.MOCSimulator(tm, results_obj,friction)

#%%
tm = tsnet.network.TransientModel(inp_file)

# Set wavespeed
tm.set_wavespeed(1200.) # m/s
# Set time options
tf = 60   # simulation period [s]
tm.set_time(tf,dt)

# Set valve closure
tc = 0.6 # valve closure period [s]
ts = 0 # valve closure start time [s]
se = 0 # end open percentage [s]
m = 1 # closure constant [dimensionless]
valve_op = [tc,ts,se,m]
tm.valve_closure('VALVE',valve_op)

# Initialize steady state simulation
t0 = 0. # initialize the simulation at 0 [s]
engine = 'DD' # demand driven simulator
tm = tsnet.simulation.Initializer(tm, t0, engine)

# Transient simulation
results_obj = 'Tnet1' # name of the object for saving simulation results
friction = 'unsteady'
tm3 = tsnet.simulation.MOCSimulator(tm, results_obj,friction)

#%%
# report results
import matplotlib.pyplot as plt
node = 'N2'
head1 = tm1.get_node(node).head
t1 = tm1.simulation_timestamps
head2 = tm2.get_node(node).head
t2 = tm2.simulation_timestamps
head3 = tm3.get_node(node).head
t3 = tm3.simulation_timestamps
fig = plt.figure(figsize=(8,5), dpi=80, facecolor='w', edgecolor='k')
plt.plot(t1, head1, 'k',label='steady', linewidth=2.5)
plt.plot(t2, head2, 'b', label='quasi-steady',  linewidth=2.5)
plt.plot(t3, head3, 'r',label='unsteady', linewidth=2.5)
plt.xlim([t1[0],t1[-1]])
plt.xlabel("Time [s]")
plt.ylabel("Pressure Head [m]")
plt.legend(loc='best')
plt.show()
fig.savefig('tnet1_unsteady_friction.pdf', format='pdf',dpi=500)
# %%
