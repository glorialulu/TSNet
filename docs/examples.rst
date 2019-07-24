====================
Example Applications
====================

Example 1 - End-valve Closure
-----------------------------

This example shows how to simulate the closure of
valve located at the boundary of a network, as shown
in :numref:`tnet1` . There are
five steps that the application woule need to take:

.. _tnet1:
.. figure:: figures/Tnet1.PNG
   :scale: 100 %
   :alt: tnet1


1. Import wmoc packae and read the Epanet .inp file.

.. literalinclude:: ../examples/Tnet1_valve_closure.py
    :lines: 1-4

2. Set wavespeed for the pipes and time options.

.. literalinclude:: ../examples/Tnet1_valve_closure.py
    :lines: 6-11

3. Set valve operation rules, including how long it takes
    to close the valve (:math: `t_c`), when to start close the
    valve (:math: `t_s`), the open percentage when the closure
    is completed (:math: `se`), and the shape of the closure
    operation curve (:math: `m`, :math: `1` stands for linear closure,
    :math: `2` stands for quadratic closure).

.. literalinclude:: ../examples/Tnet1_valve_closure.py
    :lines: 14-20

4. Compute steady state results to establish the initial
    condition for transient simulation.

.. literalinclude:: ../examples/Tnet1_valve_closure.py
    :lines: 22-25

5. Run transient simultion.

.. literalinclude:: ../examples/Tnet1_valve_closure.py
    :lines: 27-28

After the transient simulation, the results at nodes and links
will be returned and stored in the transient model (tm) instance.
The time history of head at N3 throughout the simulation can be retrived by:

.. docstring::
    >>> print(pipe.start_node_flowrate)

To plot the head results at N3:

.. literalinclude:: ../examples/Tnet1_valve_closure.py
    :lines: 31-42

yields :numref:`tnet1_node`:

.. _tnet1_node:
.. figure:: figures/tnet1_node.pdf
   :scale: 100 %
   :alt: tnet1_node

Similarily, to plot the flowrate results in pipe P2:

.. literalinclude:: ../examples/Tnet1_valve_closure.py
    :lines: 46-57

yields :numref:`tnet1_node`:

.. _tnet1_pipe:
.. figure:: figures/tnet1_pipe.pdf
   :scale: 100 %
   :alt: tnet1_pipe









