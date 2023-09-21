import tsnet
#%%
# Open an example network and create a transient model
inp_file = '/Users/luxing/Code/TSNet/examples/networks/Tnet1.inp'
tm = tsnet.network.TransientModel(inp_file)

# Set wavespeed
tm.set_wavespeed(1200.) # m/s
# Set time options
dt = 0.1  # time step [s], if not given, use the maximum allowed dt
tf = 20   # simulation period [s]
tm.set_time(tf,dt)


# Initialize steady state simulation
t0 = 0. # initialize the simulation at 0 [s]
engine = 'DD' # demand driven simulator
tm = tsnet.simulation.Initializer(tm, t0, engine)

# Add demand pulse
tc = 1 # total demand period [s]
ts = 1 # demand pulse start time [s]
tp = 0.2 # demand pulse increase time [s]
dp = 1 # demand pulse increase multiples [s]
demand_pulse = [tc,ts,tp,dp]
tm.add_demand_pulse('N2',demand_pulse)

# Transient simulation
results_obj = 'Tnet1' # name of the object for saving simulation results
tm = tsnet.simulation.MOCSimulator(tm, results_obj)
node = 'N2'
node = tm.get_node(node)
head1 = node.head

#%%
tm = tsnet.network.TransientModel(inp_file)

# Set wavespeed
tm.set_wavespeed(1200.) # m/s
# Set time options
dt = 0.1  # time step [s], if not given, use the maximum allowed dt
tf = 20   # simulation period [s]
tm.set_time(tf,dt)

# Initialize steady state simulation
t0 = 0. # initialize the simulation at 0 [s]
engine = 'DD' # demand driven simulator
tm = tsnet.simulation.Initializer(tm, t0, engine)

# Add demand pulse
tc = 1 # total demand period [s]
ts = 1 # demand pulse start time [s]
tp = 0.2 # demand pulse increase time [s]
dp = 1 # demand pulse increase multiples [s]
demand_pulse = [tc,ts,tp,dp]
tm.add_demand_pulse('N2',demand_pulse)

tc = 1 # total demand period [s]
ts = 2 # demand pulse start time [s]
tp = 0.2 # demand pulse increase time [s]
dp = 1 # demand pulse increase multiples [s]
demand_pulse = [tc,ts,tp,dp]
tm.add_demand_pulse('N4',demand_pulse)

# Transient simulation
results_obj = 'Tnet1' # name of the object for saving simulation results
tm = tsnet.simulation.MOCSimulator(tm, results_obj)
node = 'N2'
node = tm.get_node(node)
head2 = node.head


# report results
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(8,5), dpi=80, facecolor='w', edgecolor='k')
plt.plot(tm.simulation_timestamps,head1, 'k', label ='dpa =1', linewidth=2.5)
plt.plot(tm.simulation_timestamps,head2, 'r', label ='dpa =2', linewidth=2.5)
plt.xlim([tm.simulation_timestamps[0],tm.simulation_timestamps[-1]])
plt.xlabel("Time [s]")
plt.ylabel("Pressure Head [m]")
plt.legend(loc='best')
plt.show()
fig.savefig('demand_pulse_N2.pdf', format='pdf',dpi=500)
