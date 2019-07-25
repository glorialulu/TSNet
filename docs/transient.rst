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
The transient flow is govenred by the mass and monmentum conservation
equation [WYSS93]_:

.. math::
    \frac{\partial H}{\partial t} + \frac{a^2}{g} \frac{\partial V}{\partial t} - gV\sin(\alpha) = 0

    \frac{\partial V}{\partial t} + g\frac{\partial H}{\partial t} + h_f = 0

where
:math:`H` is the piezometric head,
:math:`V` is the flow velocity in the pipe,
:math:`t` is time,
:math:`a` is the wave speed,
:math:`g` is the gravity acceleration,
:math:`\alpha` is the angle from horizontal,
and :math:`h_f` represents the head loss.

Method of Characteristics (MOC)
-------------------------------
The Method of Characteristics (MOC) method is used to solve the system of
governing equations above. The essence of MOC is to transform the set of
partial differential equations to an equivalent set of ordinary differential
equations, as shown below [LAJW99]_:

.. math::
    \frac{dV}{dt} + \frac{g}{a} \frac{dH}{dt} + h_f - gV\sin(\alpha) = 0
    only when \frac{dx}{dt} = a

    \frac{dV}{dt} - \frac{g}{a} \frac{dH}{dt} + h_f - gV\sin(\alpha) = 0
    only when \frac{dx}{dt} = -a

Headloss in pipes
---------------------

Darcy-Weisbach equation.

.. math::
    h_f = \frac{f}{2gD} V |V|

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
