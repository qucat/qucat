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


class SeriesRLC(TestCaseAppended):
    """
    Series RLC circuit parameters
    """

    def parameters(self, R, L, C):
        circuit = core.Network(
            [core.C(0, 1, C), core.L(1, 2, L), core.R(0, 2, R)]
        )
        return circuit.f_k_A_chi()

    def test_frequency(self):
        C = 100e-15
        L = 10e-9
        R = 100e-9
        w, k, A, chi = self.parameters(R, L, C)
        cpx_w = (1j * C * R + np.sqrt(4 * C * L - C**2 * R**2)) / (
            2.0 * C * L
        )
        self.assertRelativelyClose(np.real(cpx_w) / 2 / np.pi, w)

    def test_frequency2(self):
        C = 1
        L = 3
        R = 0.5
        w, k, A, chi = self.parameters(R, L, C)
        cpx_w = (1j * C * R + np.sqrt(4 * C * L - C**2 * R**2)) / (
            2.0 * C * L
        )
        self.assertRelativelyClose(np.real(cpx_w) / 2 / np.pi, w)

    def test_dissipation(self):
        C = 100e-15
        L = 10e-9
        R = 1e-3
        w, k, A, chi = self.parameters(R, L, C)
        cpx_w = (1j * C * R + np.sqrt(4 * C * L - C**2 * R**2)) / (
            2.0 * C * L
        )
        self.assertRelativelyClose(2 * np.imag(cpx_w) / 2 / np.pi, k)


class Errors(TestCaseAppended):
    def revaluing_labelled_valued_component_twice(self):
        """
        Adressing last error appearing in issue #83
        """
        cir = core.Network(
            [core.L(0, 1, 1), core.C(0, 1, 1), core.R(0, 1, "R", 1)]
        )
        try:
            cir.loss_rates(R=1)
        except Exception:
            pass

        with self.assertRaises(ValueError):
            cir.loss_rates(R=1)


class Other(TestCaseAppended):
    def multiplicity_removal(self):
        circuit = self.open_gui_file("multiple_roots.txt")
        freqs = circuit.eigenfrequencies(Cd=25e-15)
        self.assertEqual(2, len(freqs))

    def test_LC_double_series_L_double_series_C(self):
        C = 1e-8
        L = 3
        circuit = core.Network(
            [
                core.C(0, 1, C * 2),
                core.C(1, 2, C * 2),
                core.L(2, 3, L / 2),
                core.L(3, 0, L / 2),
            ]
        )
        f, k, A, chi = circuit.f_k_A_chi()
        f_expected = 1 / np.sqrt(L * C) / 2 / np.pi
        self.assertRelativelyClose(f_expected, f)


class TransmonResonator(TestCaseAppended):
    def parameters(self, Cj, Lj, Cc, Cr, Lr):
        circuit = core.Network(
            [
                core.C(0, 1, Cj),
                core.J(0, 1, Lj),
                core.C(1, 2, Cc),
                core.C(0, 2, Cr),
                core.L(0, 2, Lr),
            ]
        )
        return circuit.f_k_A_chi()


class Transmon(TestCaseAppended):
    """
    Transmon circuit parameters
    """

    def parameters(self, C, Lj):
        circuit = core.Network([core.C(0, 1, C), core.J(0, 1, Lj)])
        return circuit.f_k_A_chi()

    def test_frequency(self):
        C = 100e-15
        Lj = 10e-9
        w, k, A, chi = self.parameters(C, Lj)
        self.assertRelativelyClose(1 / (np.sqrt(C * Lj) * 2.0 * pi), w)

    def test_anharmonicity(self):
        C = 100e-15
        Lj = 10e-9
        w, k, A, chi = self.parameters(C, Lj)
        self.assertRelativelyClose(e**2 / 2.0 / C / h, A[0])

    def test_phi_zpf(self):
        Cj = 100e-15
        Lj = 10e-9
        junction = core.J(0, 1, Lj)
        circuit = core.Network([core.C(0, 1, Cj), junction, core.R(0, 1, 1e6)])
        phi_0 = hbar / 2 / e
        Z = np.sqrt(Lj / Cj)
        phi_zpf = np.sqrt(hbar * Z / 2)
        self.assertRelativelyClose(
            phi_zpf / phi_0, junction.zpf(mode=0, quantity="flux")
        )

    def test_q_zpf(self):
        Cj = 100e-15
        Lj = 10e-9
        junction = core.J(0, 1, Lj)
        circuit = core.Network([core.C(0, 1, Cj), junction, core.R(0, 1, 1e6)])
        Z = np.sqrt(Lj / Cj)
        q_zpf = np.sqrt(hbar / Z / 2)
        self.assertRelativelyClose(
            q_zpf / e, np.absolute(junction.zpf(mode=0, quantity="charge"))
        )

    def test_anharmonicity_using_hamiltonian(self):
        Cj = 1e-10
        circuit = core.Network([core.C(0, 1, Cj), core.J(0, 1, 10e-9)])
        H = circuit.hamiltonian(modes=[0], taylor=4, excitations=[10])
        ee = H.eigenenergies()
        A = np.absolute((ee[1] - ee[0]) - (ee[2] - ee[1]))
        # Due to higher order terms, the mismatch with e**2/2/Cj/h is
        # (193702.3+0j) != (194712.7+0j)
        A_expected = 194712.7
        self.assertRelativelyClose(A_expected, A)

    def test_double_series_capacitor(self):
        C = 100e-15
        Lj = 10e-9
        circuit = core.Network(
            [core.C(0, 1, C * 2), core.C(1, 2, C * 2), core.J(0, 2, Lj)]
        )
        f, k, A, chi = circuit.f_k_A_chi()
        self.assertArrayRelativelyClose(
            [e**2 / 2.0 / C / h, 1 / (np.sqrt(C * Lj) * 2.0 * pi)],
            [A[0], f[0]],
        )


class SweepingParameters(TestCaseAppended):
    """
    Coupled transmon/RLC
    """

    def test_sweeping_LJ_in_fkAchi(self):
        cir = core.Network(
            [
                core.C(0, 1, 100e-15),
                core.J(0, 1, "L_J"),
                core.C(1, 2, 1e-15),
                core.C(2, 0, 100e-15),
                core.L(2, 0, 10e-9),
                core.R(2, 0, 1e6),
            ]
        )
        [cir.f_k_A_chi(L_J=x) for x in [1e-9, 2e-9]]


class TestGraphics(TestCaseAppended):
    def test_error_when_trying_to_plot_from_Network_show(self):
        circuit = core.Network([core.C(0, 1, "C"), core.J(0, 1, "Lj")])
        with self.assertRaises(TypeError):
            circuit.show()

    def test_error_when_trying_to_plot_from_Network_show_normal_modes(self):
        circuit = core.Network([core.C(0, 1, "C"), core.J(0, 1, "Lj")])
        with self.assertRaises(TypeError):
            circuit.show_normal_mode()

    def test_sweeping_CJ_array_in_zpf(self):
        C_comp = core.C(0, 1, "C_J")
        cir = core.Network(
            [
                C_comp,
                core.J(0, 1, 10e-9),
                core.C(1, 2, 1e-15),
                core.C(2, 0, 100e-15),
                core.L(2, 0, 10e-9),
                core.R(2, 0, 1e6),
            ]
        )
        self.assertRelativelyClose(
            C_comp.zpf(mode=1, quantity="charge", C_J=1.5e-9),
            C_comp.zpf(mode=1, quantity="charge", C_J=[1e-9, 1.5e-9, 3e-9])[1],
        )


class TestNetworkAnalysis(TestCaseAppended):
    def test_transfer_left_right_port_identical(self):
        """
        Trivial cases
        """
        net = core._Network([core.R(0, 1, "Z2")])
        self.assertEqual(net.transfer(0, 1, 0, 1), 1)

    def test_transfer_left_right_port_indentical_inverted(self):
        """
        Trivial cases
        """
        net = core._Network([core.R(0, 1, "Z2")])
        self.assertEqual(net.transfer(0, 1, 1, 0), -1)

    def test_transfer_voltage_divider(self):
        """
        Voltage divider, see:
        https://en.wikipedia.org/wiki/Voltage_divider
        We add an extra resistor between Vin and ground
        """

        # Compute the bridge transfer function
        net = core._Network(
            [
                core.R(0, 1, "Z2"),
                core.R(1, 2, "Z1"),
                core.R(2, 0, "Zg"),
            ]
        )
        transfer = net.transfer(0, 2, 0, 1)

        # What the tranfer function should be
        def transfer_theory(Z1, Z2, Zg):
            return Z2 / (Z1 + Z2)

        # Define some numerical values
        # for the resistors
        test_parameters = {"Z1": 1, "Z2": 2, "Zg": 3}

        self.assertRelativelyClose(
            transfer.evalf(subs=test_parameters),
            transfer_theory(**test_parameters),
        )

    def test_transfer_wheatstone(self):
        """
        Wheatstone bridge, see:
        https://en.wikipedia.org/wiki/Wheatstone_bridge
        Note that the correct definitions of the voltages V_G and V_S can be found here:
        http://www.ece.lsu.edu/ee4770/1999/lsli04.4up.pdf
        """

        # Compute the bridge transfer function
        net = core._Network(
            [
                core.R(0, 1, "R_3"),
                core.R(1, 2, "R_x"),
                core.R(2, 3, "R_2"),
                core.R(3, 0, "R_1"),
            ]
        )
        transfer = net.transfer(0, 2, 3, 1)

        # What the tranfer function should be
        def transfer_theory(R_1, R_2, R_3, R_x):
            return R_2 / (R_1 + R_2) - R_x / (R_x + R_3)

        # Define some numerical values
        # for the bridge resistors
        test_parameters = {
            "R_1": 1,
            "R_2": 2,
            "R_3": 3,
            "R_x": 4,
        }

        self.assertRelativelyClose(
            transfer.evalf(subs=test_parameters),
            transfer_theory(**test_parameters),
        )

    def test_transfer_wheatstone_all_equal(self):
        """
        Case where all resistors of the bridge are equal, and
        there is potential for these lines:
            A_lattice = (Ya + Yb)*(Yd + Yc)/(Ya*Yd-Yb*Yc)
            B_lattice = (Ya + Yb + Yc + Yd)/(Ya*Yd-Yb*Yc)
        to raise an error.
        Note: I made sure with a break-point that this code
        actually calls those two lines.
        """

        # Compute the bridge transfer function
        net = core._Network(
            [
                core.R(0, 1, "R_3"),
                core.R(1, 2, "R_x"),
                core.R(2, 3, "R_2"),
                core.R(3, 0, "R_1"),
            ]
        )
        transfer = net.transfer(0, 2, 3, 1)

        # Define some numerical values
        # for the bridge resistors
        test_parameters = {
            "R_1": 1,
            "R_2": 1,
            "R_3": 1,
            "R_x": 1,
        }

        # We just run this to check wether it raises an error
        transfer.evalf(subs=test_parameters)

    def test_transfer_square_symmetrical_lattice(self):
        """
        Square symmetrical lattice
        """

        # Compute the transfer function
        left_minus = 0
        left_plus = 1
        right_minus = 2
        right_plus = 3
        net = core._Network(
            [
                core.R(left_minus, left_plus, "Zs"),
                core.R(right_minus, right_plus, "Zl"),
                core.Admittance(left_minus, right_plus, 0),
                core.Admittance(right_minus, left_plus, 0),
                core.R(left_minus, right_minus, "Za"),
                core.R(right_plus, left_plus, "Za"),
            ]
        )
        transfer = net.transfer(left_minus, left_plus, right_minus, right_plus)

        # What the tranfer function should be
        def transfer_theory(Zs, Zl, Za):
            return Zl / (2 * Za + Zl)

        # Define some numerical values
        # for the resistors
        test_parameters = {
            "Zs": 1,
            "Zl": 2,
            "Za": 3,
        }

        self.assertRelativelyClose(
            transfer.evalf(subs=test_parameters),
            transfer_theory(**test_parameters),
        )

    def test_open_or_series_check(self):
        with self.assertRaises(ValueError):
            net = core._Network(
                [
                    core.R(0, 1, "Z"),
                    core.C(1, 2, "Z"),
                    core.J(2, 3, "Z"),
                ]
            )

    def test_shorted_LC(self):
        with self.assertRaises(ValueError):
            self.open_gui_file("shorted_circuit_3.txt")

    def test_shorted_circuit_1(self):
        with self.assertRaises(ValueError):
            self.open_gui_file("shorted_circuit_1.txt")

    def test_shorted_circuit_2(self):
        with self.assertRaises(ValueError):
            self.open_gui_file("shorted_circuit_2.txt")

    def test_connectivity_check_single_element_not_connected(self):
        with self.assertRaises(ValueError):
            net = core._Network(
                [
                    core.R(0, 1, "Z"),
                    core.C(1, 2, "Z"),
                    core.J(3, 4, "Z"),
                ]
            )

    def test_connectivity_check_subcircuit_not_connected(self):
        with self.assertRaises(ValueError):
            net = core._Network(
                [
                    core.R(0, 1, "Z"),
                    core.C(1, 2, "Z"),
                    core.J(3, 4, "Z"),
                    core.J(3, 4, "Z"),
                    core.J(4, 5, "Z"),
                ]
            )


if __name__ == "__main__":
    unittest.main()
