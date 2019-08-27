*************
About
*************

Quantum circuits constructed from Josephson junctions and superconducting electronics
are key to many quantum computing and quantum optics applications. Designing these
circuits involves calculating the Hamiltonian describing their quantum behavior. QuCAT, 
or “Quantum Circuit Analyzer Tool”, is an open-source framework to
help in this task. This open-source Python library features an intuitive graphical or 
programmatical interface to create circuits, the ability to compute their Hamiltonian, and
a set of complimentary functionalities such as calculating dissipation rates or visualizing
current flow in the circuit. QuCAT currently supports quantization in the basis of 
normal modes.

How QuCAT works
=================

For an overview of the circuit quantization method used by QuCAT and the algorithmic methods which
implement it, go to our technical paper LINK TO COME.

Authors
=======

QuCAT is currently developped and maintained by `Mario Gely <https://scholar.google.com/citations?user=Y3EbVooAAAAJ&hl=en>`_ 
in the group of `Gary Steele <http://steelelab.tudelft.nl>`_ at the University of Delft in the Netherlands.

Contact
=======

Don't hesitate to contact Mario at mario.gely@qucat.org

Contributing 
============

Your contribution is more than welcome!
You can submit pull requests `on our Github <https://github.com/mgely/qucat/>`_, or contact Mario if you want to bring big contributions to the project.

Possible extensions of the QuCAT features could include black-box impedance components to model distributed components `[1] <https://arxiv.org/abs/1204.0587>`_, 
more precisely modeling lossy circuits `[2] <https://arxiv.org/abs/1403.7341>`_, `[3] <https://arxiv.org/abs/1505.04116>`_, 
handling static offsets in flux or charge through DC sources, additional elements such as coupled inductors or 
superconducting quantum interference devices (SQUIDS) and different quantization methods, enabling for example 
quantization in the charge or flux basis. 
The latter would extend QuCAT beyond the scope of weakly-anharmonic circuits.

In terms of performance, QuCAT would benefit from delegating analytical calculations to a more efficient, 
compiled language. 

[1] arxiv.org/abs/1204.0587

[2] arxiv.org/abs/1403.7341

[3] arxiv.org/abs/1505.04116

Funding
=======

This work is supported by the European Research Council under the European Union’s H2020 program with grant agreements 681476 - QOM3D and 732894 - HOT.