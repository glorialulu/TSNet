import wmoc 

# open an example network
inp_file = 'networks/Tnet1.inp'

dt = 0.01  # time step [s]
tf = 2 # simulation peroid [s] 

# name of the valve to be close
valve_to_close = ['InlineValve']
# valve operation rule  
tc = 2 # valve closure peroid
ts = 0 # valve closure start time 
se = 0 # end open percentage 
m = 1 # closure constant
valve_op = [tc,ts,se,m]

wn, H1, V1, tt = wmoc.simulation.MOCSimulator(inp_file, dt, tf, 
                        valve_to_close, valve_op) 

# specify the pipe name to visulize results
pipe = '0'
wmoc.postprocessing.plot_head_history(pipe,H1,wn,tt)