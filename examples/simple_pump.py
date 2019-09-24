import tsnet
# Open an example network and create a transient model
inp_file = 'examples/networks/simple_pump.inp'
tm = tsnet.network.TransientModel(inp_file)

# Set wavespeed
tm.set_wavespeed(1200.) # m/s
# Set time options
dt = 0.1  # time step [s], if not given, use the maximum allowed dt
tf = 60   # simulation period [s]
tm.set_time(tf)

# Set pump shut off
tc = 1 # pump closure period
ts = 0 # pump closure start time
se = 0 # end open percentage
m = 1 # closure constant
pump_op = [tc,ts,se,m]
tm.pump_shut_off('pump', pump_op)

# Initialize steady state simulation
t0 = 0. # initialize the simulation at 0 [s]
engine = 'DD' # demand driven simulator
tm = tsnet.simulation.Initializer(tm, t0, engine)

# Transient simulation
results_obj = 'simple_pump' # name of the object for saving simulation results
tm = tsnet.simulation.MOCSimulator(tm, results_obj)

# report results
import matplotlib.pyplot as plt
node = '2'
node = tm.get_node(node)
fig1 = plt.figure(figsize=(10,4), dpi=80, facecolor='w', edgecolor='k')
plt.plot(tm.simulation_timestamps,node.head)
plt.xlim([tm.simulation_timestamps[0],tm.simulation_timestamps[-1]])
plt.title('Pressure Head at Node %s '%node)
plt.xlabel("Time [s]")
plt.ylabel("Pressure Head [m]")
plt.legend(loc='best')
plt.grid(True)
plt.show()
fig1.savefig('./docs/figures/tnet1_node.png', format='png',dpi=100)

pipe = 'p2'
pipe = tm.get_link(pipe)
fig = plt.figure(figsize=(10,4), dpi=80, facecolor='w', edgecolor='k')
plt.plot(tm.simulation_timestamps,pipe.start_node_flowrate,label='Start Node')
plt.plot(tm.simulation_timestamps,pipe.end_node_flowrate,label='End Node')
plt.xlim([tm.simulation_timestamps[0],tm.simulation_timestamps[-1]])
plt.title('Velocity of Pipe %s '%pipe)
plt.xlabel("Time [s]")
plt.ylabel("Velocity [m/s]")
plt.legend(loc='best')
plt.grid(True)
plt.show()