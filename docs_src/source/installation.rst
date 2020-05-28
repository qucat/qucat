*************
Installation
*************

.. note:: QuCAT only supports python 3

.. note:: If you don't succeed in installing QuCAT following the instructions below, ask for help on `our forum <https://groups.google.com/forum/#!forum/qucat>`_.

Installing via pip
==================

The recommended way to install qucat is via pip by opening a terminal and running

``pip install qucat``

Installing from source
======================

The latest source code is available on our Github repository

`<https://github.com/qucat/qucat>`_

To install from source, download or clone the source code, 
open a terminal and navigate to the qucat folder, and run 
``pip install .``


Requirements
============

Qucat depends on several open-source libraries. 
The following packages are currently required:

* Python 3, tested on version 3.7
* Numpy, tested on version 1.16.2
* Matplotlib, tested on version 3.0.3
* Sympy, tested on version 1.3

Generating a Hamiltonian requires

* `QuTiP <http://qutip.org/docs/latest/installation.html>`_, tested on version 4.3

Performance of Sympy and thus QuCAT is enhanced by using

* `gmpy2 <https://gmpy2.readthedocs.io/en/latest/>`_

We recommend installing python and these packages by 
downloading and installing 
`Anaconda <https://www.anaconda.com/distribution/>`_.