import wmoc
# open an example network and create a transient model
inp_file = 'examples/networks/Tnet1.inp'
tm = wmoc.network.TransientModel(inp_file)

# set wavespeed
tm.set_wavespeed(1200.) # m/s
#set time options
dt = 0.1  # time step [s], if not given, use the maximum allowed dt
tf = 60   # simulation period [s]
tm.set_time(tf,dt)

# set valve closing
tc = 0.6 # valve closure period [s]
ts = 0 # valve closure start time [s]
se = 0 # end open percentage [s]
m = 1 # closure constant [dimensionless]
valve_op = [tc,ts,se,m]
tm.valve_closure('VALVE',valve_op)

# Initialize
t0 = 0. # initialize the simulation at 0 [s]
engine = 'DD' # demand driven simulator
tm = wmoc.simulation.Initializer(tm, t0, engine)

# Transient simulation
tm = wmoc.simulation.MOCSimulator(tm)

# report results
import matplotlib.pyplot as plt

node = 'N3'
node = tm.get_node(node)
fig1 = plt.figure(figsize=(10,4), dpi=80, facecolor='w', edgecolor='k')
plt.plot(tm.simulation_timestamps,node.head)
plt.xlim([tm.simulation_timestamps[0],tm.simulation_timestamps[-1]])
plt.title('Pressure Head at Node %s '%node)
plt.xlabel("Time")
plt.ylabel("Pressure Head (m)")
plt.legend(loc='best')
plt.grid(True)
plt.show()
fig1.savefig('./docs/figures/tnet1_node.png', format='png',dpi=100)

pipe = 'P2'
pipe = tm.get_link(pipe)
fig2 = plt.figure(figsize=(10,4), dpi=80, facecolor='w', edgecolor='k')
plt.plot(tm.simulation_timestamps,pipe.start_node_flowrate,
        label='Start Node')
plt.plot(tm.simulation_timestamps,pipe.end_node_flowrate,
        label='End Node')
plt.xlim([tm.simulation_timestamps[0],tm.simulation_timestamps[-1]])
plt.title('Flowrate of Pipe %s '%pipe)
plt.xlabel("Time")
plt.ylabel("Flow rate (m^3/s)")
plt.legend(loc='best')
plt.grid(True)
plt.show()
fig2.savefig('./docs/figures/tnet1_pipe.png', format='png',dpi=100)