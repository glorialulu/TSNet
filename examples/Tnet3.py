import wmoc

# open an example network and creat a transient model
inp_file = 'examples/networks/Tnet3.inp'
tm = wmoc.network.TransientModel(inp_file)

# set wavespeed
tm.set_wavespeed(1200.)

#set time step
dt = 0.05  # time step [s], if not given, use the maximum allowed dt
tf = 20 # simulation peroid [s]
tm.set_time(tf)


# add burst
ts = 1 # burst start time
tc = 1 # time for burst to fully develop
final_burst_coeff = 0.01 # final burst coeff [ m^3/s/(m H20)^(1/2)]
tm.add_burst('JUNCTION-20', ts, tc, final_burst_coeff)

# Initialize
t0 = 0. # initialize the simulation at 0s
engine = 'WNTR' # or Epanet
tm = wmoc.simulation.Initializer(tm, t0, engine)

# Transient simulation
tm = wmoc.simulation.MOCSimulator(tm)

# report results
import matplotlib.pyplot as plt
node = 'JUNCTION-22'
node = tm.get_node(node)
fig = plt.figure(figsize=(10,4), dpi=80, facecolor='w', edgecolor='k')
plt.plot(tm.simulation_timestamps,node.head)
plt.xlim([tm.simulation_timestamps[0],tm.simulation_timestamps[-1]])
plt.title('Pressure Head at Node %s '%node)
plt.xlabel("Time")
plt.ylabel("Pressure Head (m)")
plt.legend(loc='best')
plt.grid(True)
plt.show()

pipe = 'LINK-40'
pipe = tm.get_link(pipe)
fig = plt.figure(figsize=(10,4), dpi=80, facecolor='w', edgecolor='k')
plt.plot(tm.simulation_timestamps,pipe.start_node_velocity,label='Start Node')
plt.plot(tm.simulation_timestamps,pipe.end_node_velocity,label='End Node')
plt.xlim([tm.simulation_timestamps[0],tm.simulation_timestamps[-1]])
plt.title('Velocity of Pipe %s '%pipe)
plt.xlabel("Time")
plt.ylabel("Velocity (m/s)")
plt.legend(loc='best')
plt.grid(True)
plt.show()