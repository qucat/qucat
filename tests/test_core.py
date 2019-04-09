import unittest
import Qcircuits.core as core
from math import isclose
import numpy as np
from scipy.constants import e, pi, h, hbar

class TestCaseAppended(unittest.TestCase):
    def assertRelativelyClose(self,a,b):
        a = float('%.13e'%(np.real(a))) + 1j*float('%.13e'%(np.imag(a)))
        b = float('%.13e'%(np.real(b))) + 1j*float('%.13e'%(np.imag(b)))
        self.assertEquals(a,b)

class SymbolicOperations(TestCaseAppended):

    def test_calculating_Y(self):
        C = 1.1e-14
        L = 1.3e-8
        Lj = 1.33e-9
        w = 33e9
        circuit = core.Qcircuit_NET([
            core.C(0,1,C),
            core.L(1,2,L),
            core.J(0,2,Lj)
        ])
        self.assertRelativelyClose(
            1j*(-1/Lj/w + 1/(-L*w + 1/(C*w))),
            circuit.Y_lambdified(w))

    def test_calculating_dY(self):
        # only works on resonance
        C = 1.1e-14
        L = 1.3e-8
        Lj = 1.33e-9
        w = np.sqrt(1/C/(L+Lj))
        circuit = core.Qcircuit_NET([
            core.C(0,1,C),
            core.L(1,2,L),
            core.J(0,2,Lj)
        ])
        self.assertRelativelyClose(
            1j*(1/(Lj*w**2) - (-L - 1/(C*w**2))/(1/(C*w) - L*w)**2),
            circuit.dY(w))



class StandardQuantumCircuits(TestCaseAppended):

    '''
    Transmon circuit parameters
    '''

    def transmon_parameters(self,C,Lj):
        circuit = core.Qcircuit_NET([
            core.C(0,1,C),
            core.J(0,1,Lj)
        ])
        return circuit.w_k_A_chi()

    def test_transmon_frequency(self):
        C = 100e-15
        Lj = 10e-9
        w,k,A,chi = self.transmon_parameters(C,Lj)
        self.assertRelativelyClose(1/(np.sqrt(C*Lj)*2.*pi),w)

    def test_transmon_anharmonicity(self):
        C = 100e-15
        Lj = 10e-9
        w,k,A,chi = self.transmon_parameters(C,Lj)
        self.assertRelativelyClose(e**2/2./C/h,A[0])

    '''
    Shunted Josephson ring
    '''
    
    def shunted_josephson_ring_parameters(self,C,L):
        circuit = core.Qcircuit_NET([
            core.C(0,2,C),
            core.C(1,3,C),
            core.J(0,1,L),
            core.J(1,2,L),
            core.J(2,3,L),
            core.J(3,0,L)
        ])
        return circuit.w_k_A_chi()

    def test_shunted_josephson_ring_number_of_modes_nHz(self):
        C = 2.e9
        L = 3.e11
        w,k,A,chi = self.shunted_josephson_ring_parameters(C,L)
        self.assertEqual(len(w),2,msg = f"f_res = {w}")
    def test_shunted_josephson_ring_number_of_modes_Hz(self):
        C = 2.
        L = 3. 
        w,k,A,chi = self.shunted_josephson_ring_parameters(C,L)
        self.assertEqual(len(w),2,msg = f"f_res = {w}")

    def test_shunted_josephson_ring_number_of_modes_GHz(self):
        C = 1e-13
        L = 1e-8
        w,k,A,chi = self.shunted_josephson_ring_parameters(C,L)
        self.assertEqual(len(w),2,msg = f"f_res = {w}")


class TestTesting(TestCaseAppended):
    def test_test(self):
        self.assertEqual(0,0)
        
class TestNetworkAnalysis(TestCaseAppended):

    def test_transfer_left_right_port_identical(self):
        '''
        Trivial cases
        '''
        net = core.Network([])
        self.assertEqual(net.transfer(0,1,0,1),1)

    def test_transfer_left_right_port_indentical_inverted(self):
        '''
        Trivial cases
        '''
        net = core.Network([])
        self.assertEqual(net.transfer(0,1,1,0),-1)

    def test_transfer_voltage_divider(self):
        '''
        Voltage divider, see:
        https://en.wikipedia.org/wiki/Voltage_divider
        We add an extra resistor between Vin and ground
        '''

        # Compute the bridge transfer function
        net = core.Network([
            core.R(0,1,'Z2'),
            core.R(1,2,'Z1'),
            core.R(2,0,'Zg'),
            ])
        transfer = net.transfer(0,2,0,1)

        # What the tranfer function should be
        def transfer_theory(Z1,Z2,Zg):
            return Z2/(Z1+Z2)

        # Define some numerical values
        # for the resistors
        test_parameters = {
            'Z1':1,
            'Z2':2,
            'Zg':3
        }
        
        self.assertRelativelyClose(
            transfer.evalf(subs = test_parameters),
            transfer_theory(**test_parameters))

    def test_transfer_wheatstone(self):
        '''
        Wheatstone bridge, see:
        https://en.wikipedia.org/wiki/Wheatstone_bridge
        Note that the correct definitions of the voltages V_G and V_S can be found here:
        http://www.ece.lsu.edu/ee4770/1999/lsli04.4up.pdf
        '''

        # Compute the bridge transfer function
        net = core.Network([
            core.R(0,1,'R_3'),
            core.R(1,2,'R_x'),
            core.R(2,3,'R_2'),
            core.R(3,0,'R_1'),
            ])
        transfer = net.transfer(0,2,3,1)

        # What the tranfer function should be
        def transfer_theory(R_1,R_2,R_3,R_x):
            return R_2/(R_1+R_2)-R_x/(R_x+R_3)

        # Define some numerical values
        # for the bridge resistors
        test_parameters = {
            'R_1':1,
            'R_2':2,
            'R_3':3,
            'R_x':4,
        }

        self.assertRelativelyClose(
            transfer.evalf(subs = test_parameters),
            transfer_theory(**test_parameters))

    def test_transfer_wheatstone_all_equal(self):
        '''
        Case where all resistors of the bridge are equal, and
        there is potential for these lines:
            A_lattice = (Ya + Yb)*(Yd + Yc)/(Ya*Yd-Yb*Yc)
            B_lattice = (Ya + Yb + Yc + Yd)/(Ya*Yd-Yb*Yc)
        to raise an error.
        Note: I made sure with a break-point that this code
        actually calls those two lines.
        '''

        # Compute the bridge transfer function
        net = core.Network([
            core.R(0,1,'R_3'),
            core.R(1,2,'R_x'),
            core.R(2,3,'R_2'),
            core.R(3,0,'R_1'),
            ])
        transfer = net.transfer(0,2,3,1)

        # Define some numerical values
        # for the bridge resistors
        test_parameters = {
            'R_1':1,
            'R_2':1,
            'R_3':1,
            'R_x':1,
        }

        # We just run this to check wether it raises an error
        transfer.evalf(subs = test_parameters)

    def test_transfer_square_symmetrical_lattice(self):
        '''
        Square symmetrical lattice
        '''

        # Compute the transfer function
        left_minus = 0 
        left_plus = 1
        right_minus = 2
        right_plus = 3
        net = core.Network([
            core.R(left_minus,left_plus,'Zs'),
            core.R(right_minus,right_plus,'Zl'),
            core.Admittance(left_minus,right_plus,0),
            core.Admittance(right_minus,left_plus,0),
            core.R(left_minus,right_minus,'Za'),
            core.R(right_plus,left_plus,'Za')
            ])
        transfer = net.transfer(left_minus,left_plus,right_minus,right_plus)

        # What the tranfer function should be
        def transfer_theory(Zs,Zl,Za):
            return Zl/(2*Za+Zl)

        # Define some numerical values
        # for the resistors
        test_parameters = {
            'Zs':1,
            'Zl':2,
            'Za':3,
        }

        self.assertRelativelyClose(
            transfer.evalf(subs = test_parameters),
            transfer_theory(**test_parameters))

if __name__ == '__main__':
    pass