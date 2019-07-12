WMOC
====


.. image:: https://img.shields.io/pypi/v/wmoc.svg
        :target: https://pypi.python.org/pypi/wmoc

.. image:: https://img.shields.io/travis/glorialulu/wmoc.svg
        :target: https://travis-ci.org/glorialulu/wmoc

.. image:: https://readthedocs.org/projects/wmoc/badge/?version=latest
        :target: https://wmoc.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status




WMOC conducts transient simulation using MOC method for water distribution systems.


* Free software: MIT license
* Github: https://github.com/glorialulu/wmoc.git
* Documentation: https://wmoc.readthedocs.io.


Features
--------

The WMOC is a Pyton package designed to perform transient simulation in water
distribution networks. The software includes capability to:

* Generate water network models based on .inp files 
* Operate elements, including pumps and valves, to generate transient events
* Add disruptive events including pipe bursts and leakages
* Perform transient simulation using MOC method
* Visulize results

For more information, go to https://wmoc.readthedocs.io.


Version
-------

WMOC is a ongoing research project, the current version is 0.1.0, which is 
still a pre-release. 


Cite WMOC
---------

To cite WMOC, use one of the following references:


Dependencies 
------------

WMOC is tested against Python versions 2.7, 3.4, 3.5 and 3.6. Further
using a Python distribution is recommended as they already contain (or easily
support installation of) many Python packages (e.g. SciPy, NumPy, pip, matplotlib,
etc.) that are used in the TEASER code. Two examples of those distributions are:

1. https://winpython.github.io/ WinPython comes along with a lot of Python
packages (e.g. SciPy, NumPy, pip, matplotlib, etc.)..
2. http://conda.pydata.org/miniconda.html Conda is an open source package
management  system and environment management system for installing multiple
versions of software  packages and their dependencies and switching easily
between them.

In addition, WMOC requires some specific Python packages:

1. wntr: Water Network Tool for Resilience 
  install on a python-enabled command line with `pip install wntr`
2. networkx: Network creation and manupulation engine
  install on a python-enabled command line with `pip install networkx`
3. pytest: Unit Tests engine
  install on a python-enabled command line with `pip install -U pytest`


Simple example
---------------

Simulate valve closure in a simple network

.. code:: python

# open an example network
inp_file = 'LoopedNet_pump_valve.inp'

# name of the valve to be close
valve_to_close = ['11']
# valve operation rule  
tc = 2 # valve closure peroid
ts = 0 # valve closure start time 
se = 0 # end open percentage 
m = 1 # closure constant
valve_op = [tc,ts,se,m]

wn, H1, V1, tt = wmoc.simulation.MOCSimulator(inp_file, dt, tf, 
                           valve_to_close, valve_op) 

# specify the pipe name to visulize results
pipe = '0'
wmoc.postprocessing.plot_head_history(pipe,H1,wn,tt)

License
-------

WMOC is released under the Revised BSD license. See the LICENSE.txt file.