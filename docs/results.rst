====================
Simulation Results
====================

Storage Structure
-----------------

Simulation results are returned and stored in the
:class:`tsnet.network.model.TransientModel` object for each node and link.

Node results include the following attributes:

    - Piezometric Head [m]
    - Emitter Discharge (including leakage and burst) [:math:`m^3/s`]
    - Demand Discharge [:math:`m^3/s`]

The result for each attribute is a Numpy array, the length of
which equals the total number of simulation time steps (:math:`tn`).
The array represents the time history of the simulation results at
each time steps.

For example, the head results at node 'JUNCTION-105' can be accessed by:

.. code:: python

  head = tm.nodes['JUNCTION-105'].head

To obtain the emitter discharge results at node 'JUNCTION-22':

.. code:: python

  emitter_discharge = tm.nodes['JUNCTION-22'].emitter_discharge

To collect the demand discharge results at node 'JUNCTION-16':

.. code:: python

  demand_discharge = tm.nodes['JUNCTION-16'].demand_discharge

Additionally, the time step (in second) and the time stamps (in seconds
from the start of the simulation) are also stored in the
:class:`tsnet.network.model.TransientModel` object. They can be retrieved
by:

.. code:: python

    dt = tm.time_step
    tt = tm.simulation_timestamps

The results can then be plotted with respect to the time stamps:

.. code:: python

    import matplotlib.pyplot as plt
    plt.plot(tt ,head)

result keys for links:

    - start_node_head
    - start_node_velocity
    - start_node_flowrate
    - end_node_head
    - end_node_velocity
    - end_node_flowrate



