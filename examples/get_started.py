import wmoc

# open an example network and creat a transient model
inp_file = 'examples/networks/Tnet1.inp'
tm = wmoc.network.TransientModel(inp_file)
# set wavespeed
tm.set_wavespeed(1200.)

# add leak
# t0 = 0. # add the leak at 0s
# leak_node = '1'
# emitter_coeff = 0.1
# tm.add_leak('1', emitter_coeff, t0)

# Initialize
t0 = 0. # initialize the simulation at 0s
dt = 0.01  # time step [s]
tf = 20 # simulation peroid [s]
tm = wmoc.simulation.Initializer(tm, t0, dt, tf)

# set pump operation
# tm.pump_shut_off

# seet valve operation
# valve = tm.getlink('')
# valve.closure()

# Transient simulation
# valve operation rule
tc = 2 # valve closure peroid
ts = 0 # valve closure start time
se = 0 # end open percentage
m = 1 # closure constant
valve_op = [tc,ts,se,m]
tm = wmoc.simulation.MOCSimulator(tm,'9',valve_op)

import matplotlib.pyplot as plt
pipe = tm.get_link('0')
plt.plot(tm.simulation_timestamps,pipe.start_node_head,label='Start Node')
plt.plot(tm.simulation_timestamps,pipe.end_node_head,label='End Node')
plt.xlim([tm.simulation_timestamps[0],tm.simulation_timestamps[-1]])
plt.title('Pressure Head of Pipe %s '%pipe)
plt.xlabel("Time")
plt.ylabel("Pressure Head (m)")
plt.legend(loc='best')
plt.grid(True)
plt.show()
# """
# valve = Operations('InlineValve')

# """
# # name of the valve to be close
# valve_to_close = ['InlineValve']
# # valve operation rule
# tc = 2 # valve closure peroid
# ts = 0 # valve closure start time
# se = 0 # end open percentage
# m = 1 # closure constant
# valve_op = [tc,ts,se,m]

# wn, H1, V1, tt = wmoc.simulation.MOCSimulator(inp_file, dt, tf,
#                         valve_to_close, valve_op)

# # specify the pipe name to visulize results
# pipe = '0'
# wmoc.postprocessing.plot_head_history(pipe,H1,wn,tt)