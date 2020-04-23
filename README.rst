TSNet
=======================


.. image:: https://img.shields.io/pypi/v/tsnet.svg
        :target: https://pypi.python.org/pypi/tsnet

.. image:: https://img.shields.io/travis/glorialulu/tsnet.svg
        :target: https://travis-ci.com/glorialulu/tsnet

.. image:: https://readthedocs.org/projects/tsnet/badge/?version=latest
        :target: https://tsnet.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://pepy.tech/badge/tsnet
        :target: https://pepy.tech/project/tsnet
        :alt: PyPI - Downloads
.. image:: https://img.shields.io/github/license/glorialulu/tsnet
        :alt: GitHub license
.. image:: https://img.shields.io/github/release-date-pre/glorialulu/TSNet
        :alt: GitHub (Pre-)Release Date

TSNet performs transient simulation in water networks using Method of Characteristics (MOC).


* Free software: MIT license
* GitHub: https://github.com/glorialulu/TSNet.git
* Documentation: https://tsnet.readthedocs.io.

Overview
---------

A number of commercial software for transient simulation in water
distribution systems are available in the market; however, the use of
these software for research purposes is limited. The major restriction is
due to the fact that the programs are packed as black boxes, and the source
code is not visible, thus prohibiting any changes, including modification of
existing and implementation of new elements, in the source code.
Therefore, the authors find it imperative to develop an open source package
rendering easiness for interaction, modification, and extension.

Features
--------

TSNet is a Python package designed to perform transient simulation in water
distribution networks. The software includes capabilities to:

* Create water network models based on .inp files
* Generate transient events by operating valves and pumps
* Add disruptive events including pipe bursts and leakages
* Perform transient simulation using MOC method
* Visualize results

For more information, go to https://tsnet.readthedocs.io.


Version
-------

TSNet is a ongoing research project in the University of Texas at Austin.
The current version is 0.2.0, which is still a pre-release.

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

Lu Xing and Lina Sela. " TSNet: a Python package
for transient simulations in water distribution networks".
Submitted to Advances in Engineering Software.

License
-------

TSNet is released under the MIT license. See the LICENSE.txt file.
