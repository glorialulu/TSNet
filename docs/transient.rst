========================
Transient Simulation
========================

After the steady state calculation is completed, the Method of Characteristics
(MOC) is used for solving governing transient flow equations. A transient
simulation can be run using the following code:

.. literalinclude:: ../examples/Tnet1_valve_closure.py
    :lines: 27-27

The results will then be returned to the transient model (tm) instance.


Mass and Momentum Conservation
-------------------------------
The Method of Characteristics (MOC) method is used to solve the following
governing equations transient flow [WYSS93]_.

.. math::
    \frac{\partial H}{\partial t} + a^2 /g * \frac{\partial V}{\partial t} = 0

..math::
    \frac{\partial V}{\partial t} + g*\frac{\partial H}{\partial t} + h_f = 0

where
:math:H is the piezometric head,
:math:V is the flow velocity in the pipe,
:math:a is the wave speed,
:math:g is the gravity acceleration,
and :math:h_f represents the head loss.

Method of Characteristics (MOC)
-------------------------------


Headloss in pipes
---------------------



Pressure-driven Demand
----------------------







Inner nodal boundary
--------------------




Nodal leakage boundary
----------------------




Choice of time step
-----------------------



Valve closure
--------------
Simulate valve closure



Pump Shut-off
--------------
simulate pump controlled shut-off



Burst
-----
simulate pipe burst
