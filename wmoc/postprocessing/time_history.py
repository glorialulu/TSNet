"""
The wmoc.postprocessing.time_history module contains functions 
to plot the time history of head and velocity at the start and 
end point of a pipe 
"""
from __future__ import print_function
import matplotlib.pyplot as plt

def plot_head_history(pipe,H,wn,tt):
    """Plot Head history on the start and end node of a pipe
    
    Parameters
    ----------
    pipe : str
        Name of the pipe where you want to report the head
    H : list 
        Head results 
    wn : wntr.network.model.WaterNetworkModel
        Network
    tt : list
        Simulation timestamps
    """

    pipeid = wn.links[pipe].id-1
    plt.figure(figsize=(10,4), dpi=80, facecolor='w', edgecolor='k') 
    plt.plot(tt,H[pipeid][0,:], 'b-',label='Start Node') 
    plt.plot(tt,H[pipeid][-1,:], 'r-',label='End Node') 
    plt.xlim([tt[0],tt[-1]])
    plt.title('Pressure Head of Pipe %s '%pipe)
    plt.xlabel("Time")
    plt.ylabel("Pressure Head (m)") 
    plt.legend(loc='best')
    plt.grid(True)
    plt.show()

def plot_velocity_history(pipe,V,wn,tt):
    """Plot Velocity history on the start and end node of a pipe
    
    Parameters
    ----------
    pipe : str
        Name of the pipe where you want to report the head
    V : list 
        velocity results 
    wn : wntr.network.model.WaterNetworkModel
        Network
    tt : list
        Simulation timestamps
    """

    pipeid = wn.links[pipe].id-1
    plt.figure(figsize=(10,4), dpi=80, facecolor='w', edgecolor='k')
    plt.plot(tt,V[pipeid][0,:], 'b-',label='Start Node') 
    plt.plot(tt,V[pipeid][-1,:], 'r-',label='End Node') 
    plt.xlim([tt[0],tt[-1]])
    plt.title('Velocity Head of Pipe %s ' %pipe)
    plt.xlabel("Time")
    plt.ylabel("Velocity (m/s)") 
    plt.legend(loc='best')
    plt.grid(True)
    plt.show()