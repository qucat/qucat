import sys, os

sys.path.append(os.path.join(os.path.dirname(__file__), "src"))

from core import Network, GUI, J, L, C, R


c = Network([C(0, 1, 100e-15), L(1, 2, 10e-9), R(0, 2, 0.1),])
c.Q_min = 1
print(c.eigenfrequencies())
