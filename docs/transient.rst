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
equations that apply along specific lines, i.e., characteristics lines,
as shown below [LAJW99]_:

.. math::
    \frac{dV}{dt} + \frac{g}{a} \frac{dH}{dt} + h_f - gV\sin(\alpha) = 0
    \text{only when} \frac{dx}{dt} = a

    \frac{dV}{dt} - \frac{g}{a} \frac{dH}{dt} + h_f - gV\sin(\alpha) = 0
    \text{only when} \frac{dx}{dt} = -a

The explicit MOC technique is then adopted to solve the above systems of
equations along the characteristics lines [LAJW99]_.

Headloss in Pipes
---------------------

WMOC adopts Darcy-Weisbach equation to compute head loss, regardless of the
friction method defined in the EPANet .inp file. This package computes
Darcy-Weisbach coefficients (:math:`f`) based on the head loss (:math:`{h_l}_0`)
and flow velocity (:math:`V_0`) in initial condition, using the following equation:

.. math::
    f = \frac{{h_l}_0}{(L/D)(V_0^2/2g)}

where
:math:`L` is the pipe length,
:math:`D` is the pipe diameter,
and :math:`g` is gravity acceleration.

Pressure-driven Demand
----------------------

During the transient simulation in WMOC, the demands are treated as pressure-
dependent discharge, thus indicating that the actual demands are not
equivalent to the demands defined in the .inp file.

The actual demands (:math:`D_{actual}`) are modeled based on the
instantaneous pressure, and the demand discharge coefficients,
using the following equation:

.. math::
    D_{actual} = k \sqrt{H}

where :math:`H` is the head and :math:`k` is the demand discharge coefficient,
which is calculated from the initial demand (:math:`D_0`) and head (:math:`H_0`):

.. math::
    k = D_0/ \sqrt{H_0}

It should be noted that if the head is negative, the demand flow will be
treated zero, assuming that a backflow preventer exists on each node.


Choice of Time Step
-----------------------

The determination of time step in MOC is not a trivial task. There are two
conflicting constraints that have to be satisfied simultaneously::

1.  the time step has to be the same for any pipe in the network, no matter
    how long or short they are;

2.  The Courant;s criterion has to be satisfied for each pipe as well:

.. math:: 









Valve Operation (Closure and Opening)
-------------------------------------
Simulate valve closure



Pump Operation (Shut-off and Start-up)
--------------------------------------
simulate pump controlled shut-off



Leakage
--------



Burst
-----
simulate pipe burst
