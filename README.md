## Bbox: Black box quantization in Python

Bbox is an implementation of the Black-box circuit quantization method described in https://arxiv.org/abs/1204.0587.
We provide an easy way to construct a circuit, to plot it, and to return the resonance frequencies and corresponding anharmonicities.
An example usage can be found in the jupyter notebook "Examples.ipynb"

## Requirements:

- Commonplace scientific libraries (sympy,numpy, matplotlib, ..)
- Lcapy. The best way to obtain it is to download (or clone) the source code from github: https://github.com/mph-/lcapy, and then run "python setup.py install" from the downloaded folder.