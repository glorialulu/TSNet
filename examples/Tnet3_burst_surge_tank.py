import tsnet
# open an example network and create a transient model
inp_file = 'networks/Tnet3.inp'
tm = tsnet.network.TransientModel(inp_file)

# Set wavespeed
import numpy as np
wavespeed = np.random.normal(1200., 100., size=tm.num_pipes)
tm.set_wavespeed(wavespeed)
# Set time step
tf = 20 # simulation period [s]
tm.set_time(tf)

# Add burst
ts = 1 # burst start time
tc = 1 # time for burst to fully develop
final_burst_coeff = 0.01 # final burst coeff [ m^3/s/(m H20)^(1/2)]
tm.add_burst('JUNCTION-73', ts, tc, final_burst_coeff)

# Initialize steady state simulation
t0 = 0. # initialize the simulation at 0s
engine = 'DD' # or Epanet
tm = tsnet.simulation.Initializer(tm, t0, engine)

# Transient simulation
result_obj = 'Tnet3' # name of the object for saving simulation results
tm1 = tsnet.simulation.MOCSimulator(tm,result_obj)
#%%
tm = tsnet.network.TransientModel(inp_file)

# Set wavespeed
import numpy as np
wavespeed = np.random.normal(1200., 100., size=tm.num_pipes)
tm.set_wavespeed(wavespeed)
# Set time step
tf = 20 # simulation period [s]
tm.set_time(tf)

# Add burst
ts = 1 # burst start time
tc = 1 # time for burst to fully develop
final_burst_coeff = 0.01 # final burst coeff [ m^3/s/(m H20)^(1/2)]
tm.add_burst('JUNCTION-73', ts, tc, final_burst_coeff)

# Initialize steady state simulation
t0 = 0. # initialize the simulation at 0s
engine = 'DD' # or Epanet
tm = tsnet.simulation.Initializer(tm, t0, engine)

tank_height = 10  # tank height [m]
water_height = 5  # initial water level [m]
tank_node = 'JUNCTION-30'
tank_area = 10   # tank cross sectional area [m^2]
tm.add_surge_tank(tank_node, [tank_area,tank_height,water_height], 'closed')

# Transient simulation
result_obj = 'Tnet3' # name of the object for saving simulation results
tm2 = tsnet.simulation.MOCSimulator(tm,result_obj)
#%%
# report results
import matplotlib.pyplot as plt
node1 = 'JUNCTION-16'
node2 = 'JUNCTION-20'
node3 = 'JUNCTION-30'
node4 = 'JUNCTION-45'
node5 = 'JUNCTION-90'
fig = plt.figure(figsize=(8,5), dpi=80, facecolor='w', edgecolor='k')
plt.plot(tm1.simulation_timestamps,tm1.get_node(node1).head-tm1.get_node(node1).head[0],'C0--')
#plt.plot(tm1.simulation_timestamps,tm1.get_node(node2).head-tm1.get_node(node2).head[0],'C1--')
#plt.plot(tm1.simulation_timestamps,tm1.get_node(node3).head-tm1.get_node(node3).head[0],'C2--')
plt.plot(tm1.simulation_timestamps,tm1.get_node(node4).head-tm1.get_node(node4).head[0],'C3--')
#plt.plot(tm1.simulation_timestamps,tm1.get_node(node5).head-tm1.get_node(node5).head[0],'C4--')

plt.plot(tm2.simulation_timestamps,tm2.get_node(node1).head-tm2.get_node(node1).head[0],'C0',label='JUNCTION-16')
#plt.plot(tm2.simulation_timestamps,tm2.get_node(node2).head-tm2.get_node(node2).head[0],'C1',label='JUNCTION-20')
#plt.plot(tm2.simulation_timestamps,tm2.get_node(node3).head-tm2.get_node(node3).head[0],'C2',label='JUNCTION-30')
plt.plot(tm2.simulation_timestamps,tm2.get_node(node4).head-tm2.get_node(node4).head[0],'C3',label='JUNCTION-45')
#plt.plot(tm2.simulation_timestamps,tm2.get_node(node5).head-tm2.get_node(node5).head[0],'C4',label='JUNCTION-90')


plt.xlim([tm1.simulation_timestamps[0],tm1.simulation_timestamps[-1]])
plt.title('Head on Junctions')
plt.xlabel("Time [s]")
plt.ylabel("Head [m]")
plt.legend(loc='best')
plt.show()
# fig.savefig('./docs/figures/tnet3_multi.png', format='png',dpi=100)