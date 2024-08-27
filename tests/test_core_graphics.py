import sys
import os

sys.path.append(
    os.path.join(
        os.path.join(os.path.dirname(os.path.dirname(__file__)), "src")
    )
)
import unittest
import core
from math import isclose
import numpy as np
from _constants import e, pi, h, hbar
from utils import TestCaseAppended

# Run plt.ion() to avoid hanging on plt.show() calls
import matplotlib.pyplot as plt

plt.ion()


class TestGraphics(TestCaseAppended):
    def test_error_when_trying_to_plot_from_Network_show(self):
        circuit = core.Network([core.C(0, 1, "C"), core.J(0, 1, "Lj")])
        with self.assertRaises(TypeError):
            circuit.show()

    def test_error_when_trying_to_plot_from_Network_show_normal_modes(self):
        circuit = core.Network([core.C(0, 1, "C"), core.J(0, 1, "Lj")])
        with self.assertRaises(TypeError):
            circuit.show_normal_mode()

    def test_generate_graphics(self):
        import _generate_graphics

    def test_show_transmon_RLC(self):
        cir = self.open_gui_file(
            "show_normal_mode_transmon_RLC_Lj_as_parameter.txt"
        )
        cir.show()

    def test_show_normal_mode_transmon_RLC_Lj_as_parameter(self):
        cir = self.open_gui_file(
            "show_normal_mode_transmon_RLC_Lj_as_parameter.txt"
        )
        for quantity in ["flux", "voltage", "charge", "current"]:
            cir.show_normal_mode(0, quantity, Lj=1e-9)


if __name__ == "__main__":
    unittest.main()
