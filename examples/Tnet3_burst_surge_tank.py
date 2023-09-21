import tsnet
# open an example network and create a transient model
inp_file = 'networks/Tnet3.inp'
tm = tsnet.network.TransientModel(inp_file)

# Set wavespeed
import numpy as np
wavespeed = 1200
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

# add air chamber
tank_node = 'JUNCTION-89'
tank_area = 10   # tank cross sectional area [m^2]
tank_height = 10  # tank height [m]
water_height = 5  # initial water level [m]
tm.add_surge_tank(tank_node, [tank_area,tank_height,water_height], 'closed')

# Transient simulation
result_obj = 'Tnet3' # name of the object for saving simulation results
tm2 = tsnet.simulation.MOCSimulator(tm,result_obj)

#%%
tm = tsnet.network.TransientModel(inp_file)

# Set wavespeed
import numpy as np
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

# add air chamber
tank_node = 'JUNCTION-89'
tank_area = 10   # tank cross sectional area [m^2]

tm.add_surge_tank(tank_node, [tank_area], 'open')

# Transient simulation
result_obj = 'Tnet3' # name of the object for saving simulation results
tm3 = tsnet.simulation.MOCSimulator(tm,result_obj)

#%%
# report results
import matplotlib.pyplot as plt
node1 = 'JUNCTION-16'
node2 = 'JUNCTION-20'
node3 = 'JUNCTION-30'
node4 = 'JUNCTION-45'
node5 = 'JUNCTION-90'
fig, axs = plt.subplots(1,3,figsize=(15,5), dpi=80, facecolor='w', edgecolor='k')
axs[0].plot(tm1.simulation_timestamps,tm1.get_node(node1).head-tm1.get_node(node1).head[0],'C0-',label='JUNCTION-16', linewidth=2.5)
axs[0].plot(tm1.simulation_timestamps,tm1.get_node(node2).head-tm1.get_node(node2).head[0],'C1-',label='JUNCTION-20', linewidth=2.5)
axs[0].plot(tm1.simulation_timestamps,tm1.get_node(node3).head-tm1.get_node(node3).head[0],'C2-', label='JUNCTION-30',linewidth=2.5)
axs[0].plot(tm1.simulation_timestamps,tm1.get_node(node4).head-tm1.get_node(node4).head[0],'C3-', label='JUNCTION-45',linewidth=2.5)
axs[0].plot(tm1.simulation_timestamps,tm1.get_node(node5).head-tm1.get_node(node5).head[0],'C4-',label='JUNCTION-90', linewidth=2.5)
axs[0].set_xlim([tm1.simulation_timestamps[0],tm1.simulation_timestamps[-1]])
axs[0].set_ylim([-45,15])
axs[0].set_xlabel("Time [s]")
axs[0].set_ylabel("Head change [m]")
axs[0].legend(loc='best')
axs[0].set_title('(a)')
#axs[0].show()
#fig.savefig('tnet3_wo_surge_tank.pdf', format='pdf',dpi=100)

axs[1].plot(tm2.simulation_timestamps,tm2.get_node(node1).head-tm2.get_node(node1).head[0],'C0',label='JUNCTION-16', linewidth=2.5)
axs[1].plot(tm2.simulation_timestamps,tm2.get_node(node2).head-tm2.get_node(node2).head[0],'C1',label='JUNCTION-20', linewidth=2.5)
axs[1].plot(tm2.simulation_timestamps,tm2.get_node(node3).head-tm2.get_node(node3).head[0],'C2',label='JUNCTION-30', linewidth=2.5)
axs[1].plot(tm2.simulation_timestamps,tm2.get_node(node4).head-tm2.get_node(node4).head[0],'C3',label='JUNCTION-45', linewidth=2.5)
axs[1].plot(tm2.simulation_timestamps,tm2.get_node(node5).head-tm2.get_node(node5).head[0],'C4',label='JUNCTION-90', linewidth=2.5)
axs[1].set_xlim([tm1.simulation_timestamps[0],tm1.simulation_timestamps[-1]])
axs[1].set_ylim([-45,15])
axs[1].set_xlabel("Time [s]")
#axs[1].set_ylabel("Head [m]")
axs[1].legend(loc='lower right')
axs[1].set_title('(b)')

axs[2].plot(tm3.simulation_timestamps,tm3.get_node(node1).head-tm3.get_node(node1).head[0],'C0',label='JUNCTION-16', linewidth=2.5)
axs[2].plot(tm3.simulation_timestamps,tm3.get_node(node2).head-tm3.get_node(node2).head[0],'C1',label='JUNCTION-20', linewidth=2.5)
axs[2].plot(tm3.simulation_timestamps,tm3.get_node(node3).head-tm3.get_node(node3).head[0],'C2',label='JUNCTION-30', linewidth=2.5)
axs[2].plot(tm3.simulation_timestamps,tm3.get_node(node4).head-tm3.get_node(node4).head[0],'C3',label='JUNCTION-45', linewidth=2.5)
axs[2].plot(tm3.simulation_timestamps,tm3.get_node(node5).head-tm3.get_node(node5).head[0],'C4',label='JUNCTION-90', linewidth=2.5)
axs[2].set_xlim([tm1.simulation_timestamps[0],tm1.simulation_timestamps[-1]])
axs[2].set_xlabel("Time [s]")
axs[2].set_ylim([-45,15])
#axs[1].set_ylabel("Head [m]")
axs[2].legend(loc='lower right')
axs[2].set_title('(c)')

plt.show()
plt.tight_layout()
fig.savefig('tnet3_compare_surge_tank.pdf', format='pdf', bbox_inches = 'tight',dpi=100)
