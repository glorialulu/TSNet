valve_open_period=0.5


from matplotlib import pyplot as plt
import numpy
import wntr
import tsnet


def m3s_to_gpm(m3s):
    return m3s*(60*264.172)

def m_to_psi(m):
    return m*1.4219702063247


# Restore model for each pump
wn1 = wntr.network.WaterNetworkModel('pumpcurve_1.txt')
wn2 = wntr.network.WaterNetworkModel('pumpcurve_2.txt')

#%%
                                   
# Access pump curve 1
pts_1 = wn1.get_link('supply').get_pump_curve().points
coeffs_1 = wn1.get_link('supply').get_head_curve_coefficients()
q_1 = numpy.linspace(pts_1[0][0],pts_1[-1][0])
H_1 = coeffs_1[0] - coeffs_1[1]*(q_1**coeffs_1[2])

# Access pump curve 2 
pts_2 = wn2.get_link('supply').get_pump_curve().points
coeffs_2 = wn2.get_link('supply').get_head_curve_coefficients()
q_2 = numpy.linspace(pts_2[0][0],pts_2[-1][0])
H_2 = coeffs_2[0] - coeffs_2[1]*(q_2**coeffs_2[2])

# Plot pump curves
plt.plot(m3s_to_gpm(q_1), m_to_psi(H_1), label='pumpcurve_1')
plt.plot(m3s_to_gpm(q_2), m_to_psi(H_2), label='pumpcurve_2')
plt.xlabel('gpm')
plt.ylabel('PSI')
plt.draw()
plt.legend(loc='best')

# Define valve curve
valve_curve = [(100.0, 0.0008099981681764881), (90.0, 0.00040499908408824403),
               (80.0, 0.00020249954204412202), (70.0, 0.00010124977102206101),
               (60.0, 5.394587800055411e-05), (50.0, 2.7539937718000598e-05),
               (40.0, 1.6199963363529762e-05), (30.0, 9.007179630122548e-06),
               (20.0, 5.0705885327848155e-06), (10.0, 2.7053938817094702e-06),
               (0.0, 0.0)]

# Simulate using each pump curve
time_index = []
pressure = []
flowrate = []
for input_name in ['pumpcurve_1', 'pumpcurve_2']:
    # Construct transient model
    tm = tsnet.network.TransientModel(input_name+'.txt')

    # Set wavespeed
    tm.set_wavespeed(1200) # m/s

    # Set time options
    tf = 4   # simulation period [s]
    dt = tsnet.network.discretize.max_time_step(tm)
    tm.set_time(tf, dt)

    # Set valve opening
    tc = valve_open_period # valve opening period [s]
    ts = 1 # valve opening start time [s]
    se = 1 # end open fraction (ranges from 0 to 1)
    m = 1 # open constant [dimensionless]
    tm.valve_opening('valve_3',[tc,ts,se,m], curve=valve_curve)

    # Initialize steady state simulation
    tm = tsnet.simulation.Initializer(tm, 0, 'DD')

    # Transient simulation
    tm = tsnet.simulation.MOCSimulator(tm, 'temp.res', 'steady')

    # Save results
    time_index.append(numpy.asarray(tm.simulation_timestamps, dtype=float))
    pressure.append(m_to_psi(numpy.asarray(tm.get_node('ref').head, dtype=float)))
    flowrate.append(m3s_to_gpm(numpy.asarray(tm.get_link('to_ref').start_node_flowrate, dtype=float)))

# Prep plotting
fig = plt.figure(figsize=(8,5), dpi=80, facecolor='w', edgecolor='k')

# Plot pressure
plt.clf()
plt.plot(time_index[0], pressure[0], label='pump_curve_1', linewidth=1.5)
plt.plot(time_index[1], pressure[1], label='pump_curve_2', linewidth=1.5)
plt.ylabel("Pressure [PSI]")
plt.xlim([0,4])
plt.ylim([0,100])
plt.xlabel("Time [s]")
plt.legend(loc='best')
plt.draw()
plt.savefig('pressure.png')


# Plot flow-rates
plt.clf()
plt.plot(time_index[0], flowrate[0], label='pump_curve_1', linewidth=1.5)
plt.plot(time_index[1], flowrate[1], label='pump_curve_2', linewidth=1.5)
plt.ylabel("flow-rate [gpm]")
plt.xlim([0,4])
plt.ylim([0,3])
plt.xlabel("Time [s]")
plt.legend(loc='best')
plt.draw()
plt.savefig('flowrate.png')