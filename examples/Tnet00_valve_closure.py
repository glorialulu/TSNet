import tsnet
# Open an example network and create a transient model
inp_file = 'networks/Tnet00.inp'
tm = tsnet.network.TransientModel(inp_file)

# Set wavespeed
tm.set_wavespeed(1200.) # m/s
# Set time options
dt = 0.1
tf = 60   # simulation period [s]
tm.set_time(tf,dt)

# Set valve closure
tc = 0 # valve closure period [s]
ts = 2 # valve closure start time [s]
se = 0 # end open percentage [s]
m = 1 # closure constant [dimensionless]
valve_op = [tc,ts,se,m]
tm.valve_closure('3',valve_op)

# Initialize steady state simulation
t0 = 0. # initialize the simulation at 0 [s]
engine = 'DD' # demand driven simulator
tm = tsnet.simulation.Initializer(tm, t0, engine)

# Transient simulation
results_obj = 'Tnet0' # name of the object for saving simulation results
tm1 = tsnet.simulation.MOCSimulator(tm, results_obj)