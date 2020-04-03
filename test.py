import sys, os

sys.path.append(os.path.join(os.path.dirname(__file__), "src"))

from core import Network, GUI, J, L, C, R
import numpy as np

Ci = 100e-15
Li = 10e-9
Ri = 1e6
Z0 = 2e6


def S11_ana(f):
    w = f * 2 * np.pi
    w0 = 1 / np.sqrt(Li * Ci)
    ki = Z0 / (Ri + Z0) * (1 / Ri / Ci + 1 / Z0 / Ci)
    Zin = Ri / (1 + 1j * w0 / ki * (w / w0 - w0 / w))
    return (Zin - Z0) / (Zin + Z0)


c = Network([C(0, 1, Ci), R(0, 1, Ri), L(0, 1, Li), R(0, 1, Z0, "p"),])
f, k, _, _ = c.f_k_A_chi(pretty_print=True)
freqs = np.linspace(f[0] - 3 * k[0], f[0] + 3 * k[0], 31)
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

import matplotlib.pyplot as plt

plt.plot(np.real(S11), np.imag(S11), label="qucat")
plt.plot(np.real(S11_ana(freqs)), np.imag(S11_ana(freqs)), label="analytical")
plt.legend()
plt.show()
