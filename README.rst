====
WMOC
====


.. image:: https://img.shields.io/pypi/v/WMOC.svg
        :target: https://pypi.python.org/pypi/WMOC

.. image:: https://img.shields.io/travis/glorialulu/WMOC.svg
        :target: https://travis-ci.org/glorialulu/WMOC

.. image:: https://readthedocs.org/projects/WMOC/badge/?version=latest
        :target: https://WMOC.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status




WMOC conducts transient simulation using MOC method for water distribution systems.


* Free software: MIT license
* Documentation: https://WMOC.readthedocs.io.


Features
--------

To use (with caution), simply do:
	>>> import wmoc  

	>>> wn, H1, V1, tt = WMOC.wmoc(inp_file, dt, tf, valve_to_close, pump_to_operate,
                            valve_op, pump_op, pressure_zone_bc, T,
                            leak_loc, leak_coef, burst_loc, burst_coef, burst_time, 
                            block_loc, block_A) 
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
