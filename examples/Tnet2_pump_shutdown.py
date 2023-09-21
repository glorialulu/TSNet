import tsnet
# open an example network and create a transient model
inp_file = '/Users/luxing/Code/TSNet/examples/networks/Tnet2.inp'
tm = tsnet.network.TransientModel(inp_file)

# Set wavespeed
tm.set_wavespeed(1200.)
# Set time step
tf = 20 # simulation period [s]
tm.set_time(tf)

# Set pump shut off
tc = 1 # pump closure period
ts = 1 # pump closure start time
se = 0 # end open percentage
m = 1 # closure constant
pump_op = [tc,ts,se,m]
tm.pump_shut_off('PUMP2', pump_op)

# Initialize steady state simulation
t0 = 0. # initialize the simulation at 0s
engine = 'DD' # or PPD
tm = tsnet.simulation.Initializer(tm, t0, engine)

# Transient simulation
results_obj = 'Tnet2' # name of the object for saving simulation results.head
tm = tsnet.simulation.MOCSimulator(tm,results_obj)
#%%
# report results
import matplotlib.pyplot as plt
node = 'JUNCTION-105'
node = tm.get_node(node)
fig = plt.figure(figsize=(10,4), dpi=80, facecolor='w', edgecolor='w')
plt.plot(tm.simulation_timestamps, node._head, 'k', lw=3)
plt.xlim([tm.simulation_timestamps[0],tm.simulation_timestamps[-1]])
# plt.title('Pressure Head at Node %s '%node)
plt.xlabel("Time [s]", fontsize=14)
plt.ylabel("Pressure Head [m]", fontsize=14)


# plt.legend(loc='best')
# plt.grid(True)
plt.show()
# fig.savefig('./docs/figures/tnet2_node.png', format='png',dpi=100)

pipe = 'PIPE-109'
pipe = tm.get_link(pipe)
fig = plt.figure(figsize=(10,4), dpi=80, facecolor='w', edgecolor='w')
plt.plot(tm.simulation_timestamps,pipe.start_node_velocity,label='Start Node')
plt.plot(tm.simulation_timestamps,pipe.end_node_velocity,label='End Node')
plt.xlim([tm.simulation_timestamps[0],tm.simulation_timestamps[-1]])
plt.title('Velocity of Pipe %s '%pipe)
plt.xlabel("Time [s]")
plt.ylabel("Velocity [m/s]")
plt.legend(loc='best')
plt.grid(True)
plt.show()
# fig.savefig('./docs/figures/tnet2_pipe.png', format='png',dpi=100)

