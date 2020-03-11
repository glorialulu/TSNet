import tsnet
# open an example network and create a transient model
inp_file = 'networks/Tnet2.inp'
tm = tsnet.network.TransientModel(inp_file)

# Set wavespeed
tm.set_wavespeed(1200.)
# Set time step
tf = 20 # simulation period [s]
tm.set_time(tf)

# Set pump shut off
tc = 1 # pump closure period
ts = 0 # pump closure start time
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
tm1 = tsnet.simulation.MOCSimulator(tm,results_obj)


#%%
tm = tsnet.network.TransientModel(inp_file)

# Set wavespeed
tm.set_wavespeed(1200.)
# Set time step
tf = 20 # simulation period [s]
tm.set_time(tf)

# Set pump shut off
tc = 1 # pump closure period
ts = 0 # pump closure start time
se = 0 # end open percentage
m = 1 # closure constant
pump_op = [tc,ts,se,m]
tm.pump_shut_off('PUMP2', pump_op)

# Initialize steady state simulation
t0 = 0. # initialize the simulation at 0s
engine = 'DD' # or PPD

tm = tsnet.simulation.Initializer(tm, t0, engine)
tank_node = 'JUNCTION-105'
tank_height = 10  # tank height [m]
water_height = 5  # initial water level [m]
tank_area = 10   # tank cross sectional area [m^2]
tm.add_surge_tank(tank_node, [tank_area,tank_height,water_height], 'closed')

# Transient simulation
results_obj = 'Tnet2' # name of the object for saving simulation results.head

tm2 = tsnet.simulation.MOCSimulator(tm,results_obj)

#%%

# report results
import matplotlib.pyplot as plt
node = '169'
head1 = tm1.get_node(node).head
t1 = tm1.simulation_timestamps
head2 = tm2.get_node(node).head
t2 = tm2.simulation_timestamps
fig = plt.figure(figsize=(8,5), dpi=80, facecolor='w', edgecolor='k')
plt.plot(t1, head1, 'k',label='wo surge tank',  linewidth=2.5)
plt.plot(t2, head2, 'r', label='w surge tank', linewidth=2.5)
plt.xlim([t1[0],t1[-1]])
plt.xlabel("Time [s]")
plt.ylabel("Pressure Head [m]")
plt.legend(loc='best')
plt.show()
fig.savefig('tnet2_surge_tank.pdf', format='pdf',dpi=500)
# fig.savefig('./docs/figures/tnet2_node.png', format='png',dpi=100)


