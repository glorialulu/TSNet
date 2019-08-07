import tsnet
# open an example network and create a transient model
inp_file = 'examples/networks/Tnet3.inp'
tm = tsnet.network.TransientModel(inp_file)

# set wavespeed
import numpy as np
wavespeed = np.random.normal(1200., 100., size=tm.num_pipes)
tm.set_wavespeed(wavespeed)
#set time step
tf = 20 # simulation period [s]
tm.set_time(tf)

# add leak
emitter_coeff = 0.01 # [ m^3/s/(m H20)^(1/2)]
tm.add_leak('JUNCTION-22', emitter_coeff)

# add burst
ts = 1 # burst start time
tc = 1 # time for burst to fully develop
final_burst_coeff = 0.01 # final burst coeff [ m^3/s/(m H20)^(1/2)]
tm.add_burst('JUNCTION-20', ts, tc, final_burst_coeff)

# Initialize
t0 = 0. # initialize the simulation at 0s
engine = 'DD' # or Epanet
tm = tsnet.simulation.Initializer(tm, t0, engine)

# Transient simulation
tm = tsnet.simulation.MOCSimulator(tm)

# report results
import matplotlib.pyplot as plt
node = 'JUNCTION-22'
node = tm.get_node(node)
fig = plt.figure(figsize=(10,4), dpi=80, facecolor='w', edgecolor='k')
plt.plot(tm.simulation_timestamps,node.emitter_discharge)
plt.xlim([tm.simulation_timestamps[0],tm.simulation_timestamps[-1]])
plt.title('Leak discharge at Node %s '%node)
plt.xlabel("Time [s]")
plt.ylabel("Leak discharge [m^3/s]")
plt.legend(loc='best')
plt.grid(True)
plt.show()
fig.savefig('./docs/figures/tnet3_leak.png', format='png',dpi=100)

node = 'JUNCTION-20'
node = tm.get_node(node)
fig = plt.figure(figsize=(10,4), dpi=80, facecolor='w', edgecolor='k')
plt.plot(tm.simulation_timestamps,node.emitter_discharge)
plt.xlim([tm.simulation_timestamps[0],tm.simulation_timestamps[-1]])
plt.title('Burst discharge at Node %s '%node)
plt.xlabel("Time [s]")
plt.ylabel("Burst discharge [m^3/s]")
plt.legend(loc='best')
plt.grid(True)
plt.show()
fig.savefig('./docs/figures/tnet3_burst.png', format='png',dpi=100)


pipe = 'LINK-40'
pipe = tm.get_link(pipe)
fig = plt.figure(figsize=(10,4), dpi=80, facecolor='w', edgecolor='k')
plt.plot(tm.simulation_timestamps,pipe.start_node_velocity,label='Start Node')
plt.plot(tm.simulation_timestamps,pipe.end_node_velocity,label='End Node')
plt.xlim([tm.simulation_timestamps[0],tm.simulation_timestamps[-1]])
plt.title('Velocity of Pipe %s '%pipe)
plt.xlabel("Time [s]")
plt.ylabel("Velocity [m/s]")
plt.legend(loc='best')
plt.grid(True)
plt.show()
fig.savefig('./docs/figures/tnet3_pipe.png', format='png',dpi=100)

node1 = tm.get_node('JUNCTION-8')
node2 = tm.get_node('JUNCTION-16')
node3 = tm.get_node('JUNCTION-45')
node4 = tm.get_node('JUNCTION-90')
fig = plt.figure(figsize=(10,4), dpi=80, facecolor='w', edgecolor='k')
plt.plot(tm.simulation_timestamps,node1.head,label='JUNCTION-8')
plt.plot(tm.simulation_timestamps,node2.head,label='JUNCTION-16')
plt.plot(tm.simulation_timestamps,node3.head,label='JUNCTION-45')
plt.plot(tm.simulation_timestamps,node4.head,label='JUNCTION-90')
plt.xlim([tm.simulation_timestamps[0],tm.simulation_timestamps[-1]])
plt.title('Head on Junctions')
plt.xlabel("Time [s]")
plt.ylabel("Head [m]")
plt.legend(loc='best')
plt.grid(True)
plt.show()
fig.savefig('./docs/figures/tnet3_multi.png', format='png',dpi=100)