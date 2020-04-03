import sys, os

sys.path.append(os.path.join(os.path.dirname(__file__), "src"))

from core import Network, GUI, J, L, C, R
import numpy as np

c = Network([C(0, 1, 100e-15), L(0, 1, 10e-9), C(1, 2, 10e-15), R(0, 2, 50, "p"),])
f, k, _, _ = c.f_k_A_chi(pretty_print=True)
freqs = np.linspace(f[0] - 3 * k[0], f[0] + 3 * k[0], 30)
S11 = []
for fd in freqs:
    S11.append(
        c.S(
            "p",
            "p",
            fd,
            -180,
            drive_phase=0,
            power_unit="dBm",
            temperature=0,
            temperature_unit="kelvin",
            modes="all",
            taylor=2,
            excitations=4,
        )
    )
print(S11)

import matplotlib.pyplot as plt

plt.plot(freqs, np.real(S11))
plt.plot(freqs, np.imag(S11))
plt.show()
