"""
The wmoc.network.geometry read in the geometry defined by EPANet
.inp file, and assign aditional parameters needed in transient
simution later in wmoc.

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
from wmoc.network.discretize import discretization, max_time_step
from wmoc.network.control import (
    valveclosing,
    valveopening,
    pumpclosing,
    pumpopening,
    burstsetting
)

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
        self.simulation_peroid = 0.
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
            node.emitter_coeff = 0.
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

    def set_wavespeed(self, wavespeed=1200.):
        """Set wave speed for pipes in the network

        Parameters
        ----------
        wavespeed : float or int or list, optional
            If given as float or int, set the value as wavespeed
            for all pipe; If given as list set the corresponding
            value to each pipe, by default 1200.
        """

        if isinstance(wavespeed,float):
            # if wavespeed is a float, assign it to all pipes
            wavev = wavespeed * np.ones((self.num_pipes, 1))
        elif isinstance(wavespeed, (list,tuple,np.ndarray)):
            # if wavespeed is a list, assign each elements
            # to the respective pipes.
            if len(wavespeed) == self.num_pipes:
                wavev = wavespeed
            else:
                raise ValueError('The length of the wavespeed \
                input does not equal number of pipes. ')
        else:
            raise ValueError('Wavespeed should be a float or a list')

        # assign wave speed to each pipes
        i= 0
        for _, pipe in self.pipes():
            pipe.wavev = wavev[i]
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
        self.simulation_peroid = tf
        self = discretization(self, dt)
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
        burst_node.burst_coeff = burstsetting(self.time_step, self.simulation_peroid,
                                                ts, tc, final_burst_coeff)
        burst_node.burst_status = True


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
        valve.operation_rule = valveclosing(self.time_step, self.simulation_peroid, rule)

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
        valve.operation_rule = valveopening(self.time_step, self.simulation_peroid, rule)

    def pump_shut_off(self, name, rule):
        """Set pump shut off rule

        Parameters
        ----------
        name : str
            The name of the pump to shut off
        rule : list
            Contains paramtes to defie valve operation rule
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
        pump.operation_rule = pumpclosing(self.time_step, self.simulation_peroid, rule)

    def pump_start_up(self, name, rule):
        """Set pump start up rule

        Parameters
        ----------
        name : str
            The name of the pump to shut off
        rule : list
            Contains paramtes to defie valve operation rule
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
        pump.operation_rule = pumpopening(self.time_step, self.simulation_peroid, rule)








