import wmoc

# open an example network and creat a transient model
inp_file = 'examples/networks/Tnet2_v2.inp'
tm = wmoc.network.TransientModel(inp_file)

# set wavespeed
tm.set_wavespeed(1200.)

#set time step
dt = 0.05  # time step [s], if not given, use the maximum allowed dt
tf = 20 # simulation peroid [s]
tm.set_time(tf)

# add leak
# leak_node = '2'
# emitter_coeff = 0.1 #[ m^3/s/(m H20)^(1/2)]
# tm.add_leak(leak_node, emitter_coeff)

# set pump shut off
# tc = 2 # pump closure peroid
# ts = 0 # pump closure start time
# se = 0.001 # end open percentage
# m = 1 # closure constant
# pump_op = [tc,ts,se,m]
# tm.pump_shut_off('335', pump_op)

# set pump start up
tc = 2 # valve opening peroid
ts = 0 # valve opening start time
se = 1 # end open percentage
m = 1 # closure constant
pump_op = [tc,ts,se,m]
tm.pump_start_up('335', pump_op)



# Initialize
t0 = 0. # initialize the simulation at 0s
engine = 'WNTR' # or Epanet
tm = wmoc.simulation.Initializer(tm, t0, engine)

# Transient simulation
tm = wmoc.simulation.MOCSimulator(tm)

# report results
import matplotlib.pyplot as plt

node = '61'
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

pipe = '329'
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