====================
Example Applications
====================

Example 1 - End-valve Closure
-----------------------------

This example shows how to simulate the closure of
valve located at the boundary of a network. The first example
network is shown below in :numref:`tnet1`. It comprises 9 pipes,
5 junctions, one reservoir, 3 closed loops,and one valve located
at the downstream end of the system. There are five steps that the
application would need to take:

.. _tnet1:
.. figure:: figures/Tnet1.PNG
   :scale: 100 %
   :alt: tnet1


1.  Import wmoc package and read the Epanet .inp file.

.. literalinclude:: ../examples/Tnet1_valve_closure.py
    :lines: 1-4

2.  Set wave speed for all pipes to be :math: `1200m/s`,
    time step to be :math: `0.1m/s`, and simulation period
    to be :math: `60s`.

.. literalinclude:: ../examples/Tnet1_valve_closure.py
    :lines: 6-11

3.  Set valve operation rules, including how long it takes
    to close the valve (:math: `t_c`), when to start close the
    valve (:math: `t_s`), the open percentage when the closure
    is completed (:math: `se`), and the shape of the closure
    operation curve (:math: `m`, :math: `1` stands for linear closure,
    :math: `2` stands for quadratic closure).

.. literalinclude:: ../examples/Tnet1_valve_closure.py
    :lines: 13-19

4.  Compute steady state results to establish the initial
    condition for transient simulation.

.. literalinclude:: ../examples/Tnet1_valve_closure.py
    :lines: 21-24

5. Run transient simulation.

.. literalinclude:: ../examples/Tnet1_valve_closure.py
    :lines: 26-27

After the transient simulation, the results at nodes and links
will be returned and stored in the transient model (tm) instance.
The time history of head at N3 throughout the simulation can be retrieved by:

>>> print(tm.links['P2'].start_node_flowrate)

To plot the head results at N3:

.. literalinclude:: ../examples/Tnet1_valve_closure.py
    :lines: 32-43

yields :numref:`tnet1_node`:

.. _tnet1_node:
.. figure:: figures/tnet1_node.png
   :width: 600
   :alt: tnet1_node

   Head at node N3.

Similarly, to plot the flow rate results in pipe P2:

.. literalinclude:: ../examples/Tnet1_valve_closure.py
    :lines: 45-59

yields :numref:`tnet1_pipe`:

.. _tnet1_pipe:
.. figure:: figures/tnet1_pipe.png
   :width: 600
   :alt: tnet1_pipe

   Flow rate at the start and end node of pipe P2.



Example 2 - Pump operations
---------------------------

This example illustrates how the package models a controlled pump shutdown
transient event, where the pump speed is ramped down.This example
network is shown below in :numref:`tnet2`. It comprises 113 pipes,
91 junctions, one valve, two pumps, two reservoir, and three tanks.
There are five steps that the application would need to take:

.. _tnet2:
.. figure:: figures/Tnet2.PNG
   :scale: 100 %
   :alt: tnet2

1.  Import wmoc package and read the Epanet .inp file.

.. literalinclude:: ../examples/Tnet2_pump_shutdown.py
    :lines: 1-4

2.  Set wave speed for all pipes to be :math: `1200m/s` and
    simulation period to be :math: `60s`. Use suggested time
    step.

.. literalinclude:: ../examples/Tnet2_pump_shutdown.py
    :lines: 6-10

3.  Set pump operation rules, including how long it takes
    to shutdown the pump (:math: `t_c`), when to the shut-down starts
    (:math: `t_s`), the pump speed multiplier value when the shut-down
    is completed (:math: `se`), and the shape of the shut-down 
    operation curve (:math: `m`, :math: `1` stands for linear closure,
    :math: `2` stands for quadratic closure).

.. literalinclude:: ../examples/Tnet2_pump_shutdown.py
    :lines: 12-18

4.  Compute steady state results to establish the initial
    condition for transient simulation.

.. literalinclude:: ../examples/Tnet2_pump_shutdown.py
    :lines: 20-23

5. Run transient simulation.

.. literalinclude:: ../examples/Tnet2_pump_shutdown.py
    :lines: 25-26

After the transient simulation, the results at nodes and links
will be returned and stored in the transient model (tm) instance.
The actual demand discharge at JUNCTION-105 throughout the simulation
can be retrieved by:

>>> print(tm.nodes['JUNCTION-105'].demand_discharge)

To plot the head results at JUNCTION-105:

.. literalinclude:: ../examples/Tnet2_pump_shutdown.py
    :lines: 31-42

yields :numref:`tnet2_node`:

.. _tnet1_node:
.. figure:: figures/tnet2_node.png
   :width: 600
   :alt: tnet2_node

   Head at node JUNCTION-105.

Similarly, to plot the velocity results in PIPE-109:

.. literalinclude:: ../examples/Tnet1_valve_closure.py
    :lines: 44-56

yields :numref:`tnet2_pipe`:

.. _tnet2_pipe:
.. figure:: figures/tnet1_pipe.png
   :width: 600
   :alt: tnet2_pipe

   Velocity at the start and end node of PIPE-109.


Example 2 - Burst
---------------------------




