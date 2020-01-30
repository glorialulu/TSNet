"""
The tsnet.network.geometry read in the geometry defined by EPANet
.inp file, and assign additional parameters needed in transient
simulation later in tsnet.

"""

from __future__ import print_function
import wntr
from wntr.network.elements import LinkStatus
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import logging
import warnings
from wntr.network import WaterNetworkModel
from tsnet.network.discretize import (
    discretization, max_time_step,
    discretization_N, max_time_step_N)
from tsnet.network.control import (
    valveclosing,
    valveopening,
    pumpclosing,
    pumpopening,
    burstsetting
)
from tsnet.postprocessing import detect_cusum

logger = logging.getLogger(__name__)

class TransientModel (WaterNetworkModel):
    """ Transient model class.
    Parameters
    -------------------
    inp_file_name: string
        Directory and filename of EPANET inp file to load into the
        WaterNetworkModel object.
    """

    def __init__ (self, inp_file):
        super().__init__(inp_file)
        self.simulation_timestamps = []
        self.time_step = 0.
        self.simulation_period = 0.
        self.initial_velocity = []
        self.initial_head = []
        # assign ID to each links, start from 1.
        i =1
        for _, link in self.links():
            link.id = i
            i+=1

        # assign ID to each links, start from 1.
        i =1
        for _, node in self.nodes():
            node.id = i
            node.leak_status = False
            node.burst_status = False
            node.blockage_status = False
            node.emitter_coeff = 0.
            node.block_per = 0.
            i+=1     ## Graph the network

        # calculate the slope and area for each pipe
        for _, pipe in self.pipes():
            pipe.area = pipe.diameter**2. * np.pi / 4.
            try:
                theta = np.sin(np.arctan(pipe.end_node.elevation -
                    pipe.start_node.elevation)/pipe.length)
            except:
                theta = 0.0
            pipe.theta = theta

        # set operating default value as False
        for _, link in self.links():
            link.operating = False

    def set_wavespeed(self, wavespeed=1200, pipes=None):
        """Set wave speed for pipes in the network

        Parameters
        ----------
        wavespeed : float or int or list, optional
            If given as float or int, set the value as wavespeed
            for all pipe; If given as list set the corresponding
            value to each pipe, by default 1200.
        pipes : str or list, optional
            The list of pipe to define wavespeed,
            by default all pipe in the network.
        """
        generator = 0
        if pipes == None :
            generator = 1
            pipes = self.pipes()
            num_pipes = self.num_pipes
        else:
            pipes = [self.get_link(pipe) for pipe in list(pipes)]
            num_pipes = len(pipes)

        if isinstance(wavespeed,(float,int)):
            # if wavespeed is a float, assign it to all pipes
            wavespeed = wavespeed * np.ones(num_pipes)
        elif isinstance(wavespeed, (list,tuple,np.ndarray)):
            # if wavespeed is a list, assign each elements
            # to the respective pipes.
            if not len(wavespeed) == num_pipes:
                raise ValueError('The length of the wavespeed \
                input does not equal number of pipes. ')
        else:
            raise ValueError('Wavespeed should be a float or a list')

        # assign wave speed to each pipes
        i= 0
        if generator == 1:
            for _, pipe in pipes:
                pipe.wavev = wavespeed[i]
                i+=1
        else:
            for pipe in pipes:
                pipe.wavev = wavespeed[i]
                i+=1

    def set_roughness(self,roughness, pipes=None):
        """Set roughness coefficient for pipes in the network

        Parameters
        ----------
        roughness : float or int or list
            If given as float or int, set the value as roughness
            for all pipe; If given as list set the corresponding
            value to each pipe. Make sure to define it using the
            same method (H-W or D-W) as defined in .inp file.
        pipes : str or list, optional
            The list of pipe to define roughness coefficient,
            by default all pipe in the network.
        """
        generator = 0
        if pipes == None :
            generator = 1
            pipes = self.pipes()
            num_pipes = self.num_pipes
        else:
            pipes = [self.get_link(pipe) for pipe in list(pipes)]
            num_pipes = len(pipes)

        if isinstance(roughness,(float,int)):
            # if roughness is a float, assign it to all mentioned pipes
            roughness = roughness * np.ones(num_pipes)
        elif isinstance(roughness, (list,tuple,np.ndarray)):
            # if roughness is a list, assign each elements
            # to the respective pipes.
            if not len(roughness) == num_pipes:
                raise ValueError('The length of the roughness \
                input does not equal number of input pipes. ')
        else:
            raise ValueError('Roughness should be a float or a list')

        # assign roughness to each input pipes
        i= 0
        if generator == 1:
            for _, pipe in pipes:
                pipe.roughness = roughness[i]
                i+=1
        else:
            for pipe in pipes:
                pipe.roughness = roughness[i]
                i+=1

    def set_time(self, tf, dt=None):
        """Set time step and duration for the simulation.

        Parameters
        ----------
        tf : float
            Simulation period
        dt : float, optional
            time step, by default maximum allowed dt
        """
        if dt == None:
            dt = max_time_step(self)
        self.simulation_period = tf
        self = discretization(self, dt)
        print('Simulation time step %.5f s' % self.time_step)

    def set_time_N(self, tf, N=2):
        """Set time step and duration for the simulation.

        Parameters
        ----------
        tf : float
            Simulation period
        N : integer
            Number of segments in the critical pipe
        """
        dt = max_time_step_N(self,N)
        self.simulation_period = tf
        self = discretization_N(self, dt)
        print('Simulation time step %.5f s' % self.time_step)

    def add_leak(self, name, coeff):
        """Add leak to the transient model

        Parameters
        ----------
        name : str, optional
            The name of the leak nodes, by default None
        coeff : list or float, optional
            Emitter coefficient at the leak nodes, by default None
        """

        leak_node = self.get_node(name)
        leak_node.emitter_coeff += coeff
        leak_node.leak_status = True

    def add_burst(self, name, ts, tc, final_burst_coeff):
        """Add leak to the transient model

        Parameters
        ----------
        name : str
            The name of the leak nodes, by default None
        ts : float
            Burst start time
        tc : float
            Time for burst to fully develop
        final_burst_coeff : list or float
            Final emitter coefficient at the burst nodes
        """

        burst_node = self.get_node(name)
        burst_node.burst_coeff = burstsetting(self.time_step, self.simulation_period,
                                                ts, tc, final_burst_coeff)
        burst_node.burst_status = True

    def add_blockage(self, name, percentage):
        """Add blockage to the transient model

        Parameters
        ----------
        name : str
            The name of the blockage nodes, by default None
        percentage : list or float
            The percentage of the blockage flow discharge
        """
        blockage_node = self.get_node(name)
        blockage_node.block_per = percentage
        blockage_node.block_status = True

    def valve_closure(self, name, rule):
        """Set valve closure rule

        Parameters
        ----------
        name : str
            The name of the valve to close
        rule : list
            Contains paramters to define valve operation rule
            rule = [tc,ts,se,m]
            tc : the duration takes to close the valve [s]
            ts : closure start time [s]
            se : final open percentage [s]
            m  : closure constant [unitless]
        """

        valve = self.get_link(name)
        if valve.link_type.lower() != 'valve':
            raise RuntimeError('The name of valve to operate is not associated with a vale')

        if valve.status.name == 'Closed':
            warnings.warn("Valve %s is already closed in its initial setting. \
The initial setting has been changed to open to perform the closure." %name)
            valve.status = LinkStatus.Open

        valve.operating = True
        valve.operation_rule = valveclosing(self.time_step, self.simulation_period, rule)

    def valve_opening(self, name, rule):
        """Set valve opening rule

        Parameters
        ----------
        name : str
            The name of the valve to close
        rule : list
            Contains paramters to define valve operation rule
            rule = [tc,ts,se,m]
            tc : the duration takes to open the valve [s]
            ts : opening start time [s]
            se : final open percentage [s]
            m  : closure constant [unitless]
        """
        valve = self.get_link(name)
        if valve.link_type.lower() != 'valve':
            raise RuntimeError('The name of valve to operate is not associated with a vale')

        if valve.initial_status.name == 'Open' or valve.initial_status.name == 'Active':
            warnings.warn("Valve %s is already open in its initial setting. \
The initial setting has been changed to closed to perform the opening." %name)
            valve.status = LinkStatus.Closed

        valve.operating = True
        valve.operation_rule = valveopening(self.time_step, self.simulation_period, rule)

    def pump_shut_off(self, name, rule):
        """Set pump shut off rule

        Parameters
        ----------
        name : str
            The name of the pump to shut off
        rule : list
            Contains paramaters to define valve operation rule
            rule = [tc,ts,se,m]
            tc : the duration takes to close the pump [s]
            ts : closure start time [s]
            se : final open percentage [s]
            m  : closure constant [unitless]
        """
        pump = self.get_link(name)

        if pump.link_type.lower() != 'pump':
            raise RuntimeError('The name of pump to operate is not associated with a pump')

        if pump.initial_status.name == 'Closed':
            warnings.warn("Pump %s is already closed in its initial setting. \
The initial setting has been changed to open to perform the closure." %name)
            pump.status= LinkStatus.Open
        pump.operating = True
        pump.operation_rule = pumpclosing(self.time_step, self.simulation_period, rule)

    def pump_start_up(self, name, rule):
        """Set pump start up rule

        Parameters
        ----------
        name : str
            The name of the pump to shut off
        rule : list
            Contains paramaters to define valve operation rule
            rule = [tc,ts,se,m]
            tc : the duration takes to close the valve [s]
            ts : closure start time [s]
            se : final open percentage [s]
            m  : closure constant [unitless]
        """
        pump = self.get_link(name)
        if pump.link_type.lower() != 'pump':
            raise RuntimeError('The name of pump to operate is not associated with a pump')

        # Turn the pump on and run initial calculation
        # to get the nominal flow and head
        pump.status = LinkStatus.Open
        sim = wntr.sim.EpanetSimulator(self)
        results = sim.run_sim()
        pump.nominal_flow = results.link['flowrate'].loc[0,name]
        node1 = self.links[name].start_node.name
        node2 = self.links[name].end_node.name
        pump.nominal_pump_head = abs(results.node['head'].loc[0,node1]-
                        results.node['head'].loc[0,node2])

        # Turn the pump back to closed
        pump.status = LinkStatus.Closed
        pump.operating = True
        pump.operation_rule = pumpopening(self.time_step, self.simulation_period, rule)

    def detect_pressure_change(self, name, threshold, drift, show=False, ax=None):
        """Detect pressure change in simulation results

        Parameters
        ----------
        name : str
            The name of the node
        threshold : positive number, optional (default = 1)
            amplitude threshold for the change in the data.
        drift : positive number, optional (default = 0)
            drift term that prevents any change in the absence of change.
        show : bool, optional (default = True)
            True (1) plots data in matplotlib figure, False (0) don't plot.
        ax : a matplotlib.axes.Axes instance, optional (default = None).
        """
        time = self.simulation_timestamps
        x = self.get_node(name).head
        ta, tf, amp = detect_cusum(time, x, threshold, drift,
                 show, ax=None)
        ta = [time[i] for i in ta]
        tf = [time[i] for i in tf]
        print ('%s changes detected in pressure results on node %s' %(len(ta), name))

        return ta, tf, list(amp)

    def plot_node_head(self, name, ax=None):
        """Detect pressure change in simulation results

        Parameters
        ----------
        name : str or list
            The name of node
        ax : a matplotlib.axes.Axes instance, optional (default = None).
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print('matplotlib is not available.')
        else:
            if ax is None:
                fig, ax = plt.subplots(1, 1, figsize=(8,4),dpi=100, facecolor='w', edgecolor='k')
        if not type(name) is list:
            name = [name]
        nodes = [self.get_node(i) for i in name]
        time = self.simulation_timestamps

        for i,node in enumerate(nodes):
            ax.plot(time,node.head,lw=2,label=name[i])
        plt.xlim([self.simulation_timestamps[0],self.simulation_timestamps[-1]])
        # plt.title('Pressure Head at Node(s) ')
        plt.xlabel("Time [s]", fontsize=14)
        plt.ylabel("Pressure Head [m]", fontsize=14)
        plt.legend(loc='best', framealpha=.5, numpoints=1)
        plt.grid(False)
        plt.show()








