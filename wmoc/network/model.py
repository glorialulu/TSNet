"""
The wmoc.network.geometry read in the geometry defined by EPANet
.inp file, and assign aditional parameters needed in transient
simution later in wmoc.

"""

from __future__ import print_function
import wntr
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import logging
from wntr.network import WaterNetworkModel
import collections
from wmoc.network.control import valvesetting, pumpsetting,burstsetting

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
            node.leaking = False
            node.bursting = False
            node.emitter_coeff = 0.
            i+=1     ## Graph the network

        # calculate the slope for each pipe
        for _, pipe in self.pipes():
            try:
                theta = np.sin(np.arctan(pipe.end_node.elevation -
                    pipe.start_node.elevation)/pipe.length)
            except:
                theta = 0.0
            pipe.theta = theta

        # set operating defualt value as False
        for _, link in self.links():
            link.operating = False

    def set_wavespeed(self, wavespeed=1200.):
        """Set wave speed for pipes in the network

        Parameters
        ----------
        wavespeed : float or int or list, optional
            If given as float or int, set the value as wavespeed
            for all pipe; If given as list set the correspondiing
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
        leak_node.leaking = True

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
                                                ts,tc, final_burst_coeff)
        burst_node.bursting = True


    def valve_closure(self, name, rule):
        """Set valve closure rule

        Parameters
        ----------
        name : str
            The name of the valve to close
        rule : list
            Contains paramtes to defie valve operation rule
            valve_op = [tc,ts,se,m]
            tc : the duration takes to close the valve [s]
            ts : closure start time [s]
            se : final open percentage [s]
            m  : closure constant [unitless]
        """

        valve = self.get_link(name)
        valve.operating = True
        valve.operation_rule = valvesetting(self.time_step, self.simulation_peroid, rule)

    def pump_shut_off(self, name, rule):
        """Set pump shut off rule

        Parameters
        ----------
        name : str
            The name of the pump to shut off
        rule : list
            Contains paramtes to defie valve operation rule
            valve_op = [tc,ts,se,m]
            tc : the duration takes to close the valve [s]
            ts : closure start time [s]
            se : final open percentage [s]
            m  : closure constant [unitless]
        """
        pump = self.get_link(name)
        pump.operating = True
        pump.operation_rule = pumpsetting(self.time_step, self.simulation_peroid, rule)





