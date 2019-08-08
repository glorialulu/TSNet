=====================================
Software Conventions and Limitations
=====================================

Units
------

All data in TSNet is stored in the following SI (International System) units:

* Length = :math:`m`
* Diameter = :math:`m`
* Water pressure = :math:`m`
  (this assumes a fluid density of 1000 :math:`kg/m^3`)
* Elevation = :math:`m`
* Mass = :math:`kg`
* Time = :math:`s`
* Demand = :math:`m^3/s`
* Velocity = :math:`m/s`
* Acceleration = :math:`g` (1 :math:`g` = 9.8 :math:`m/s^2`)
* Volume = :math:`m^3`

If the unit system specified in .inp file is US units,
it will be converted to SI unit in the simulation process.
When setting up analysis in TSNet, all input values
should be specified in SI units.
All simulation results are also stored in SI units.


Assumptions and Limitations
----------------------------

Current software limitations are noted:

*   Demands on the start and end nodes of pumps and valves are not supported.
    If demands are defined on these nodes in the .inp file, they will be
    ignored in transient simulation, and the simulation results may
    not be accurate due to discrepancies between the initial conditions
    and the first step in transient simulation. Warnings will be printed.

*   Multi-branch junctions on the start and end nods of pumps and valves
    are not supported. It is assumed that valves and pumps are connected
    by series pipes.

*   Demands are pressure dependent during transient simulation.

*   Friction coefficients are converted to D-W coefficients based
    on initial conditions.

*   Pipe bursts and leaks occur only on the node.







