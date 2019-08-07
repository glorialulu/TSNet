====================
Simulation Results
====================

Storage Structure
-----------------

Simulation results are returned and stored in the
:class:`tsnet.network.model.TransientModel` object for each node and link.

Node results include the following attributes:

    - Piezometric head [m]
    - Emitter discharge (including leakage and burst) [:math:`m^3/s`]
    - Demand discharge [:math:`m^3/s`]

Link results include the following attributes:

    - Head at start node [m]
    - Flow velocity at start node [:math:`m^3/s`]
    - Flowrate at start node [:math:`m^3/s`]
    - Head at end node [m]
    - Flow velocity at end node [:math:`m^3/s`]
    - Flowrate at end node [:math:`m^3/s`]


The result for each attribute is a Numpy array, representing the time
history of the simulation results, the length of
which equals the total number of simulation time steps (:math:`tn`).

For example, the results of head, emitter discharge and demand discharge
at node 'JUNCTION-105' can be accessed by:

.. code:: python

  node = tm.get_node['JUNCTION-105']
  head = node.head
  emitter_discharge = node.emitter_discharge
  demand_discharge = node.demand_discharge

To obtain the results on pipe 'LINK-40':

.. code:: python

  pipe = tm.get_link('LINK-40')
  start_head = pipe.start_node_head
  end_head = pipe.end_node_head
  start_velocity = pipe.start_node_velocity
  end_velocity = pipe.end_node_velocity
  start_flowrate = pipe.start_node_flowrate
  end_flowrate = pipe.end_node_flowrate


Time Step and Time Stamps
-------------------------

Additionally, the time step (in second) and the time stamps (in seconds
from the start of the simulation) are also stored in the
:class:`tsnet.network.model.TransientModel` object. They can be retrieved
by:

.. code:: python

    dt = tm.time_step
    tt = tm.simulation_timestamps

The results can then be plotted with respect to the time stamps using
**matplotlib** or any other preferred package, as shown in :numref:`tnet2_node`:

.. code:: python

    import matplotlib.pyplot as plt
    plt.plot(tt ,head)

.. _tnet2_node:
.. figure:: figures/tnet2_node.png
   :width: 600
   :alt: tnet2_node

Results Retrieval
------------------

The :class:`tsnet.network.model.TransientModel` object, including
the information of the network, operation rules, and the simulated results,
is saved in the file **results.obj**, located in the current folder.

To retrieve the results from a previously completed simulation, one can read
the :class:`tsnet.network.model.TransientModel` object from the
**results.obj** file and access results from the objet just like shown above:

.. code:: python

    import pickle
    file = open('results.obj', 'rb')
    tm = pickle.load(file)
