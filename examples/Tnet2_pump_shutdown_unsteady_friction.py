import tsnet
# open an example network and create a transient model
inp_file = 'networks/Tnet2.inp'
tm = tsnet.network.TransientModel(inp_file)

# Set wavespeed
tm.set_wavespeed(1200.)
# Set time step
tf = 10 # simulation period [s]
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
friction ='steady'
tm1 = tsnet.simulation.MOCSimulator(tm,results_obj,friction)


#%%
tm = tsnet.network.TransientModel(inp_file)

# Set wavespeed
tm.set_wavespeed(1200.)
# Set time step
tf = 10 # simulation period [s]
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
friction ='quasi-steady'
tm2 = tsnet.simulation.MOCSimulator(tm,results_obj,friction)

#%%
tm = tsnet.network.TransientModel(inp_file)

# Set wavespeed
tm.set_wavespeed(1200.)
# Set time step
tf = 10 # simulation period [s]
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
friction ='unsteady'
tm3 = tsnet.simulation.MOCSimulator(tm,results_obj,friction)
#%%

# report results
import matplotlib.pyplot as plt
node = 'JUNCTION-105'
head1 = tm1.get_node(node).head
t1 = tm1.simulation_timestamps
head2 = tm2.get_node(node).head
t2 = tm2.simulation_timestamps
head3 = tm3.get_node(node).head
t3 = tm3.simulation_timestamps
fig = plt.figure(figsize=(8,5), dpi=80, facecolor='w', edgecolor='k')
plt.plot(t1, head1, 'k',label='steady',  linewidth=2.5)
plt.plot(t2, head2, 'b', label='quasi-steady', linewidth=2.5)
plt.plot(t3, head3, 'r',label='unsteady',  linewidth=2.5)
plt.xlim([t1[0],t1[-1]])
plt.xlabel("Time [s]")
plt.ylabel("Pressure Head [m]")
plt.legend(loc='best')
plt.show()
fig.savefig('tnet2_unsteady_friction.pdf', format='pdf',dpi=500)
# fig.savefig('./docs/figures/tnet2_node.png', format='png',dpi=100)


