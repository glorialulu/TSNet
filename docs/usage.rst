===============
Getting Started
===============

To use tsnet in a project, open a Python console and import the package::

    import tsnet

Simple example
---------------

A simple example, Tnet1_valve_closure.py is included in the examples folder.
This example demonstrates how to:

* Import tsnet

* Generate a transient model

* Set wave speed

* Set time step and simulation period

* Perform initial condition calculation

* Define valve closure rule

* Run transient simulation and save results to .obj file

* Plot simulation results

.. literalinclude:: ../examples/Tnet1_valve_closure.py
    :lines: 2-42


Three additional EPANET INP files and example files are also included
in the TSNet examples repository in the examples folder.
Example networks range from a simple 8-node network to a 126-node network.
