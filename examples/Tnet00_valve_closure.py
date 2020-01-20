import tsnet
# Open an example network and create a transient model
inp_file = 'networks/Tnet00.inp'
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
valve_op = [tc,ts,se,m]
tm.valve_closure('3',valve_op)

# Initialize steady state simulation
t0 = 0. # initialize the simulation at 0 [s]
engine = 'DD' # demand driven simulator
tm = tsnet.simulation.Initializer(tm, t0, engine)

# Transient simulation
results_obj = 'Tnet0' # name of the object for saving simulation results
tm = tsnet.simulation.MOCSimulator(tm, results_obj)

node = '3'
node = tm.get_node(node)

dt=tm.simulation_timestamps
norm_head_nl = node.head - node.head[0]

#%%
# report results
import matplotlib.pyplot as plt
fig1 = plt.figure(figsize=(10,4), dpi=80, facecolor='w', edgecolor='k')
plt.plot(tm.simulation_timestamps,norm_head_nl,'k--')
plt.xlim([tm.simulation_timestamps[0],tm.simulation_timestamps[-1]])
plt.xlabel("Time [s]")
plt.ylabel("Pressure Head Change [m]")
#plt.legend(loc='best')
plt.grid(False)
plt.show()
# fig1.savefig('./docs/figures/tnet1_node.png', format='png',dpi=100)
