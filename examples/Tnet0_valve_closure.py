import tsnet
# Open an example network and create a transient model
inp_file = '/Users/luxing/Code/TSNet/examples/networks/Tnet0.inp'
tm = tsnet.network.TransientModel(inp_file)

# Set wavespeed
tm.set_wavespeed(1200.) # m/s
# Set time options
dt = 0.01
tf = 25   # simulation period [s]
tm.set_time(tf,dt)

# Set valve closure
tc = 0 # valve closure period [s]
ts = 0 # valve closure start time [s]
se = 0 # end open percentage [s]
m = 1 # closure constant [dimensionless]
valve_op = [tc, ts, se, m]
tm.valve_closure('3',valve_op)

# Initialize steady state simulation
t0 = 0. # initialize the simulation at 0 [s]
engine = 'PDD' # demand driven simulator
tm = tsnet.simulation.Initializer(tm, t0, engine)

# Transient simulation
results_obj = 'Tnet0' # name of the object for saving simulation results
friction = 'steady'
tm1 = tsnet.simulation.MOCSimulator(tm, results_obj, friction)

#%%
tm = tsnet.network.TransientModel(inp_file)

# Set wavespeed
tm.set_wavespeed(1200.) # m/s
# Set time options
tm.set_time(tf,dt)
tm.valve_closure('3',valve_op)

# Initialize steady state simulation
tm = tsnet.simulation.Initializer(tm, t0, engine)

# Transient simulation
friction = 'quasi-steady'
tm2 = tsnet.simulation.MOCSimulator(tm, results_obj, friction)

#%%
tm = tsnet.network.TransientModel(inp_file)

# Set wavespeed
tm.set_wavespeed(1200.) # m/s
# Set time options

tm.set_time(tf,dt)


# Set valve closure
tm.valve_closure('3',valve_op)

# Initialize steady state simulation
tm = tsnet.simulation.Initializer(tm, t0, engine)

# Transient simulation
friction = 'unsteady'
tm3 = tsnet.simulation.MOCSimulator(tm, results_obj, friction)
#%%
# report results
import matplotlib.pyplot as plt
node = '2'
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
fig.savefig('tnet0_unsteady_friction.pdf', format='pdf',dpi=500)
