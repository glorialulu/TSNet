====================
Example Applications
====================

Example 1 - End-valve Closure
-----------------------------

This example shows how to simulate the closure of
valve located at the boundary of a network, as shown
in |tnet1| . There are
three steps that the application woule need to take:

.. |tnet1| image:: figures/Tnet1.PNG

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








