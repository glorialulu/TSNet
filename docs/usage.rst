===========
Get started
===========

To use wmoc in a project, open a Python console and import the package::

    import wmoc


Simple example
---------------

A simple example, get_startted.py is included in the examples folder.
Thie example demonstrate how to:

* Import wmoc

* Gnerate a transient model

* Set wave speed

* Set time step and simulation peroid

* Perform initial condition calculation

* Define valve closure rule

* Run transient simulation

* Plot simulation results

.. literalinclude:: ../examples/Tnet1_valve_closure.py
    :lines: 2-27

Additional examples are included in the wmoc documentation.

Three EPANET INP files and example files are also included in the wmoc
repository in the examples folder. Example networks range from a simple
8 node network to a 126 node network.
