import sys
import os
sys.path.append(os.path.join(os.path.join(os.path.dirname(os.path.dirname(__file__)),'src')))
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
    '''
    Series RLC circuit parameters
    '''

    def parameters(self,R,L,C):
        circuit = core.Network([
            core.C(0,1,C),
            core.L(1,2,L),
            core.R(0,2,R)
        ])
        return circuit.f_k_A_chi()

    def test_frequency(self):
        C = 100e-15
        L = 10e-9
        R = 100e-9
        w,k,A,chi = self.parameters(R,L,C)
        cpx_w = (1j*C*R + np.sqrt(4*C*L - C**2*R**2))/(2.*C*L)
        self.assertRelativelyClose(np.real(cpx_w)/2/np.pi,w)

    def test_frequency2(self):
        C = 1
        L = 3
        R = 0.5
        w,k,A,chi = self.parameters(R,L,C)
        cpx_w = (1j*C*R + np.sqrt(4*C*L - C**2*R**2))/(2.*C*L)
        self.assertRelativelyClose(np.real(cpx_w)/2/np.pi,w)


    def test_dissipation(self):
        C = 100e-15
        L = 10e-9
        R = 100e-9
        w,k,A,chi = self.parameters(R,L,C)
        cpx_w = (1j*C*R + np.sqrt(4*C*L - C**2*R**2))/(2.*C*L)
        self.assertRelativelyClose(np.imag(cpx_w)/2/np.pi,k)

class Other(TestCaseAppended):

    def multiplicity_removal(self):
        circuit = self.open_gui_file('multiple_roots.txt')
        freqs = circuit.eigenfrequencies(Cd = 25e-15)
        self.assertEqual(2,len(freqs))

    def test_LC_double_series_L_double_series_C(self):
        C = 1e-8
        L = 3
        circuit = core.Network([
            core.C(0,1,C*2),
            core.C(1,2,C*2),
            core.L(2,3,L/2),
            core.L(3,0,L/2)
        ])
        f,k,A,chi = circuit.f_k_A_chi()
        f_expected = 1/np.sqrt(L*C)/2/np.pi
        self.assertRelativelyClose(f_expected,f)

    def test_0_value_in_list_kw(self):
        cir = core.Network([
            core.C(0,1,'C'),
            core.J(0,1,1e-9)])
        with self.assertRaises(ValueError):
            cir.f_k_A_chi(C = np.linspace(1,0,101))
            
    def test_0_value_in_single_kw(self):
        cir = core.Network([
            core.C(0,1,'C'),
            core.J(0,1,1e-9)])
        with self.assertRaises(ValueError):
            cir.f_k_A_chi(C = 0)
        
class TransmonResonator(TestCaseAppended):

    def parameters(self,Cj,Lj,Cc,Cr,Lr):
        circuit = core.Network([
            core.C(0,1,Cj),
            core.J(0,1,Lj),
            core.C(1,2,Cc),
            core.C(0,2,Cr),
            core.L(0,2,Lr)
        ])
        return circuit.f_k_A_chi()

class Transmon(TestCaseAppended):
    '''
    Transmon circuit parameters
    '''

    def parameters(self,C,Lj):
        circuit = core.Network([
            core.C(0,1,C),
            core.J(0,1,Lj)
        ])
        return circuit.f_k_A_chi()

    def test_frequency(self):
        C = 100e-15
        Lj = 10e-9
        w,k,A,chi = self.parameters(C,Lj)
        self.assertRelativelyClose(1/(np.sqrt(C*Lj)*2.*pi),w)

    def test_anharmonicity(self):
        C = 100e-15
        Lj = 10e-9
        w,k,A,chi = self.parameters(C,Lj)
        self.assertRelativelyClose(e**2/2./C/h,A[0])

    def test_phi_zpf(self):
        Cj = 100e-15
        Lj = 10e-9
        junction = core.J(0,1,Lj)
        circuit = core.Network([
            core.C(0,1,Cj),
            junction,
            core.R(0,1,1e6)
        ])
        phi_0 = hbar/2/e
        Z = np.sqrt(Lj/Cj)
        phi_zpf = np.sqrt(hbar*Z/2)
        self.assertRelativelyClose(phi_zpf/phi_0,junction.zpf(mode=0,quantity = 'flux'))

        
    def test_q_zpf(self):
        Cj = 100e-15
        Lj = 10e-9
        junction = core.J(0,1,Lj)
        circuit = core.Network([
            core.C(0,1,Cj),
            junction,
            core.R(0,1,1e6)
        ])
        Z = np.sqrt(Lj/Cj)
        q_zpf = np.sqrt(hbar/Z/2)
        self.assertRelativelyClose(q_zpf/e,np.absolute(junction.zpf(mode=0,quantity = 'charge')))

    def test_anharmonicity_using_hamiltonian(self):
        Cj = 1e-10
        circuit = core.Network([
            core.C(0,1,Cj),
            core.J(0,1,10e-9)
        ])
        H = circuit.hamiltonian(modes = [0],taylor = 4,excitations = [10])
        ee = H.eigenenergies()
        A = np.absolute((ee[1]-ee[0])-(ee[2]-ee[1]))
        # Due to higher order terms, the mismatch with e**2/2/Cj/h is
        # (193702.3+0j) != (194712.7+0j)
        A_expected = 194712.7
        self.assertRelativelyClose(A_expected,A)

    def test_double_series_capacitor(self):
        C = 100e-15
        Lj = 10e-9
        circuit = core.Network([
            core.C(0,1,C*2),
            core.C(1,2,C*2),
            core.J(0,2,Lj)
        ])
        f,k,A,chi = circuit.f_k_A_chi()
        self.assertArrayRelativelyClose([e**2/2./C/h,1/(np.sqrt(C*Lj)*2.*pi)],[A[0],f[0]])

class SweepingParameters(TestCaseAppended):
    '''
    Coupled transmon/RLC
    '''
    def test_sweeping_LJ_forloop_in_fkAchi(self):
        cir = core.Network([
            core.C(0,1,100e-15),
            core.J(0,1,'L_J'),
            core.C(1,2,1e-15),
            core.C(2,0,100e-15),
            core.L(2,0,10e-9),
            core.R(2,0,1e6)
            ])
        [cir.f_k_A_chi(L_J=x) for x in [1e-9,2e-9]]

    def test_sweeping_LJ_nparray_in_fkAchi(self):
        cir = core.Network([
            core.C(0,1,100e-15),
            core.J(0,1,'L_J'),
            core.C(1,2,1e-15),
            core.C(2,0,100e-15),
            core.L(2,0,10e-9),
            core.R(2,0,1e6)
            ])
        f,k,A,chi = cir.f_k_A_chi(L_J=np.array([1e-9,2e-9]))
        self.assertRelativelyClose(A[1,0],cir.anharmonicities(L_J=1e-9)[1])

    def test_sweeping_LJ_array_in_fkAchi(self):
        cir = core.Network([
            core.C(0,1,100e-15),
            core.J(0,1,'L_J'),
            core.C(1,2,1e-15),
            core.C(2,0,100e-15),
            core.L(2,0,10e-9),
            core.R(2,0,1e6)
            ])
        f,k,A,chi = cir.f_k_A_chi(L_J=[1e-9,2e-9])
        self.assertRelativelyClose(chi[1,0,1],cir.kerr(L_J=2e-9)[1,0])

    def test_sweeping_LJ_CJ_array_in_fkAchi(self):
        cir = core.Network([
            core.C(0,1,'C_J'),
            core.J(0,1,'L_J'),
            core.C(1,2,1e-15),
            core.C(2,0,100e-15),
            core.L(2,0,10e-9),
            core.R(2,0,1e6)
            ])
        f,k,A,chi = cir.f_k_A_chi(L_J=[1e-9,2e-9],C_J=[1e-9,2e-9])
        self.assertRelativelyClose(f[1,0],cir.eigenfrequencies(L_J=1e-9,C_J=1e-9)[1])
    def test_sweeping_LJ_CJ_array_in_fkAchi_incompatible_sizes(self):
        cir = core.Network([
            core.C(0,1,'C_J'),
            core.J(0,1,'L_J'),
            core.C(1,2,1e-15),
            core.C(2,0,100e-15),
            core.L(2,0,10e-9),
            core.R(2,0,1e6)
            ])
        with self.assertRaises(ValueError):
            cir.f_k_A_chi(L_J=[1e-9,2e-9],C_J=[1e-9,2e-9,3e-9])
    def test_sweeping_LJ_CJ_array_in_fkAchi_show(self):
        cir = core.Network([
            core.C(0,1,'C_J'),
            core.J(0,1,'L_J'),
            core.C(1,2,1e-15),
            core.C(2,0,100e-15),
            core.L(2,0,10e-9),
            core.R(2,0,1e6)
            ])
        
        with self.assertRaises(ValueError):
            cir.show(L_J=[1e-9,2e-9],C_J=[1e-9,2e-9,3e-9])

    def test_sweeping_LJ_CJ_array_in_fkAchi_show_normal_mode(self):
        cir = core.Network([
            core.C(0,1,'C_J'),
            core.J(0,1,'L_J'),
            core.C(1,2,1e-15),
            core.C(2,0,100e-15),
            core.L(2,0,10e-9),
            core.R(2,0,1e6)
            ])
        
        with self.assertRaises(ValueError):
            cir.show(L_J=[1e-9,2e-9],C_J=[1e-9,2e-9,3e-9])

    def test_sweeping_LJ_CJ_array_in_fkAchi_show_normal_mode(self):
        cir = core.Network([
            core.C(0,1,'C_J'),
            core.J(0,1,'L_J'),
            core.C(1,2,1e-15),
            core.C(2,0,100e-15),
            core.L(2,0,10e-9),
            core.R(2,0,1e6)
            ])
        
        with self.assertRaises(ValueError):
            cir.show(L_J=[1e-9,2e-9],C_J=[1e-9,2e-9,3e-9])

        
    def test_sweeping_CJ_array_in_zpf(self):
        C_comp = core.C(0,1,'C_J')
        cir = core.Network([
            C_comp,
            core.J(0,1,10e-9),
            core.C(1,2,1e-15),
            core.C(2,0,100e-15),
            core.L(2,0,10e-9),
            core.R(2,0,1e6)
            ])
        self.assertRelativelyClose(
             C_comp.zpf(mode = 1, quantity = 'charge', C_J=1.5e-9),
            C_comp.zpf(mode = 1, quantity = 'charge', C_J=[1e-9,1.5e-9,3e-9])[1])

class TestNetworkAnalysis(TestCaseAppended):


    def test_open_or_series_check(self):
        with self.assertRaises(ValueError):
            net = core._Network([
                core.R(0,1,'Z'),
                core.C(1,2,'Z'),
                core.J(2,3,'Z'),
                ])

    def test_shorted_LC(self):
        with self.assertRaises(ValueError):
            self.open_gui_file('shorted_circuit_3.txt')
    def test_shorted_circuit_1(self):
        with self.assertRaises(ValueError):
            self.open_gui_file('shorted_circuit_1.txt')
    def test_shorted_circuit_2(self):
        with self.assertRaises(ValueError):
            self.open_gui_file('shorted_circuit_2.txt')

    def test_connectivity_check_single_element_not_connected(self):
        with self.assertRaises(ValueError):
            net = core._Network([
                core.R(0,1,'Z'),
                core.C(1,2,'Z'),
                core.J(3,4,'Z'),
                ])
    def test_connectivity_check_subcircuit_not_connected(self):
        with self.assertRaises(ValueError):
            net = core._Network([
                core.R(0,1,'Z'),
                core.C(1,2,'Z'),
                core.J(3,4,'Z'),
                core.J(3,4,'Z'),
                core.J(4,5,'Z'),
                ])

if __name__ == "__main__":
    unittest.main()