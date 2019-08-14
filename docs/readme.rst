Introduction to TSNet
=======================
TSNet performs transient simulation in water networks using Method of Characteristics (MOC).


* Free software: MIT license
* Github: https://github.com/glorialulu/TSNet.git
* Documentation: https://tsnet.readthedocs.io.

Overview
---------
Hydraulic transients in water distribution networks (WDNs),
typically induced by pipe bursts, valve operations, and pump operations,
can disturb the steady-state flow conditions by introducing extreme pressure
variability and imposing abrupt internal pressure force
onto the pipeline systems [WOLB05]_.
These disturbances have been identified as one of the major contributing factors
in the many pipe deterioration and catastrophic failure in WDNs [RERS15]_,
thereby wasting a significant amount of treated water and creating unexpected
possibilities of contamination intrusion [ASCE17]_.
Consequently, transient simulation, as a prominent approach to
understand and predict the behavior of hydraulic transients,
has become an essential requirement for ensuring the distribution safety and
improving the efficiency in the process of design and operation of WDNs.
In addition to improving design and operation of WDNs,
various other transient-based applications, such as network calibration,
leak detection, sensor placement, and condition assessment,
has also enhanced the popularity and necessity of transient simulation

Acknowledgedly, a number of commercial software for transient simulation
in water distribution systems are available in the market;
however, the use of these software for research purposes is limited.
The major restriction is due to the fact that the programs are packed
as black boxes, and the source code is not visible,
thus prohibiting any changes, including modification of
existing and implementation of new elements, in the source code.
Additionally, the commercial software was designed to perform only
single transient simulations and do not have the capabilities to automate or
run multiple transient simulations.
Users are required to modify the boundary conditions using the GUI,
perform the simulation, and manually record the responses to changes
in the various conditions,
which significantly complicated the research process.

There is a clear gap that currently available simulation software
are not suitable for many research applications beyond the
conventional design purposes.
Hence, the motivation of this work is two-fold:

1.  Provide users with open source and freely available python code
and package for simulating transients in water distribution systems
that can be integrated with other case specific applications,
e.g. sensor placement and event detection; and

2.  Encourage users and developers to further develop and
extend the transient model.


Features
--------

TSNet is a Python package designed to perform transient simulation in water
distribution networks. The software includes capability to:

* Create transient models based on EPANET INP files
* Operating valves and pumps
* Add disruptive events including pipe bursts and leakages
* Perform transient simulation using Method of characteristics (MOC) techniques
* Visualize results

For more information, go to https://tsnet.readthedocs.io.


Version
-------

TSNet is a ongoing research project in the University of Texas at Austin.
The current version is 0.1.0, which is still a pre-release.

Contact
-------

* Lu Xing, the University of Texas at Austin, xinglu@utexas.edu
* Lina Sela, the University of Texas at Austin, linasela@utexas.edu

Disclaimer
----------

No warranty, expressed or implied, is made as to the correctness of the
results or the suitability of the application.


Cite TSNet
-----------

To cite TSNet, use one of the following references:


License
-------

TSNet is released under the MIT license. See the LICENSE.txt file.
