import tsnet
# Open an example network and create a transient model
tm = tsnet.network.TransientModel('networks/Tnet1.inp')

# Set wavespeed
tm.set_wavespeed(1200.) # m/s

# Set time options
tf = 20   # simulation period [s]
tm.set_time(tf)

# Set valve closure
ts = 5 # valve closure start time [s]
tc = 1 # valve closure period [s]
se = 0 # end open percentage [s]
m = 2 # closure constant [dimensionless]
tm.valve_closure('VALVE',[tc,ts,se,m])

# Initialize steady state simulation
t0=0
tm = tsnet.simulation.Initializer(tm,t0)

# Transient simulation
tm = tsnet.simulation.MOCSimulator(tm)

# report results
node = ['N2','N3']
tm.plot_node_head(node)