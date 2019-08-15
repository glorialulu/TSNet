.. highlight:: shell

============
Installation
============

Setup Python Environment
------------------------------

TSNet is tested against Python versions 2.7, 3.5 and 3.6.
It can be installed on Windows, Linux, and Mac OS X operating systems.
Python distributions, such as Anaconda, are recommended to manage the Python
environment as they already contain (or easily support installation of) many
Python packages (e.g. SciPy, NumPy, pandas, pip, matplotlib, etc.) that are
used in the TSNet package.  For more information on Python package
dependencies, see :ref:`Dependencies`.

Stable Release (for users)
--------------------------

To install TSNet, run this command in your terminal:

.. code-block:: console

    $ pip install tsnet

This is the preferred method to install tsnet, as it will always install the
most recent stable release.

If you don't have `pip`_ installed, this `Python installation guide`_ can guide
you through the process.

.. _pip: https://pip.pypa.io
.. _Python installation guide: http://docs.python-guide.org/en/latest/starting/installation/


From Sources (for developers)
-----------------------------

The sources for TSNet can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/glorialulu/tsnet

Or download the `tarball`_:

.. code-block:: console

    $ curl  -OL https://github.com/glorialulu/tsnet/tarball/master

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ python setup.py install


.. _Github repo: https://github.com/glorialulu/tsnet
.. _tarball: https://github.com/glorialulu/tsnet/tarball/master

.. _Dependencies:

Dependencies
------------

Requirements for TSNet include Python (2.7, 3.5, or 3.6) along with
several Python packages.
The following Python packages are required::

1.  Numpy [VaCV11]_: the fundamental package needed for scientific computing with Python
    included in Anaconda distribution
    http://www.numpy.org/

2.  Matplotlib [Hunt07]_: Python 2D plotting library
    included in Anaconda distribution
    http://matplotlib.org/

3.  NetworkX [HaSS08]_: Network creation and manipulation engine,
    install on a python-enabled command line with `pip install wntr`
    https://networkx.github.io/

4.  WNTR [WNTRSi]_: Water Network Tool for Resilience
    install on a python-enabled command line with `pip install wntr`
    http://wntr.readthedocs.io

5.  pytest: Unit Tests engine
    install on a python-enabled command line with `pip install -U pytest`
    https://docs.pytest.org/en/latest/

