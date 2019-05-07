import unittest
import qucat.src.core as core
from math import isclose
import numpy as np
from scipy.constants import e, pi, h, hbar

def cutoff_digits(f,digits):
    float_format = '%%.%de'%digits
    return float(float_format%(np.real(f))) + 1j*float(float_format%(np.imag(f)))

class TestCaseAppended(unittest.TestCase):

    def assertRelativelyClose(self,a,b,digits = 6):
        a = cutoff_digits(a,digits)
        b = cutoff_digits(b,digits)
        self.assertEqual(a,b)

    def assertArrayRelativelyClose(self,a,b,digits = 6):
        a = np.array(a)
        b = np.array(b)
        self.assertTrue(a.shape==b.shape,msg = f'Arrays do not have the same dimension {a.shape}!={b.shape}')
        for index,_ in np.ndenumerate(a):
            a_comp = cutoff_digits(a[index],digits)
            b_comp = cutoff_digits(b[index],digits)
            self.assertEqual(a_comp,b_comp,
                    msg = f'Components with index {index} do not match {a_comp}!={b_comp}')

class StandardQuantumCircuits(TestCaseAppended):
    '''
    Series RLC circuit parameters
    '''

    def series_RLC_parameters(self,R,L,C):
        circuit = core.Network([
            core.C(0,1,C),
            core.L(1,2,L),
            core.R(0,2,R)
        ])
        return circuit.f_k_A_chi()

    def test_series_RLC_frequency(self):
        C = 100e-15
        L = 10e-9
        R = 100e-9
        w,k,A,chi = self.series_RLC_parameters(R,L,C)
        cpx_w = (1j*C*R + np.sqrt(4*C*L - C**2*R**2))/(2.*C*L)
        self.assertRelativelyClose(np.real(cpx_w)/2/np.pi,w)

    def test_series_RLC_frequency2(self):
        C = 1
        L = 3
        R = 0.5
        w,k,A,chi = self.series_RLC_parameters(R,L,C)
        cpx_w = (1j*C*R + np.sqrt(4*C*L - C**2*R**2))/(2.*C*L)
        self.assertRelativelyClose(np.real(cpx_w)/2/np.pi,w)

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

    def test_series_RLC_dissipation(self):
        C = 100e-15
        L = 10e-9
        R = 100e-9
        w,k,A,chi = self.series_RLC_parameters(R,L,C)
        cpx_w = (1j*C*R + np.sqrt(4*C*L - C**2*R**2))/(2.*C*L)
        self.assertRelativelyClose(np.imag(cpx_w)/2/np.pi,k)

    '''
    Transmon circuit parameters
    '''

    def transmon_parameters(self,C,Lj):
        circuit = core.Network([
            core.C(0,1,C),
            core.J(0,1,Lj)
        ])
        return circuit.f_k_A_chi()

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

    def test_transmon_phi_zpf(self):
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

        
    def test_transmon_q_zpf(self):
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

    def test_transmon_anharmonicity_using_hamiltonian(self):
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

    def test_transmon_double_series_capacitor(self):
        C = 100e-15
        Lj = 10e-9
        circuit = core.Network([
            core.C(0,1,C*2),
            core.C(1,2,C*2),
            core.J(0,2,Lj)
        ])
        f,k,A,chi = circuit.f_k_A_chi()
        self.assertArrayRelativelyClose([e**2/2./C/h,1/(np.sqrt(C*Lj)*2.*pi)],[A[0],f[0]])

    '''
    Shunted Josephson ring
    '''
    
    def shunted_josephson_ring_parameters(self,C,L):
        circuit = core.Network([
            core.C(0,2,C),
            core.C(1,3,C),
            core.J(0,1,L),
            core.J(1,2,L),
            core.J(2,3,L),
            core.J(3,0,L)
        ])
        return circuit.f_k_A_chi()

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

    def test_shunted_josephson_ring_frequency_0(self):
        C = 1e-13
        L = 1e-8
        w,k,A,chi = self.shunted_josephson_ring_parameters(C,L)
        self.assertRelativelyClose(w[0],1/np.sqrt(L*C)/2./np.pi)

    def test_shunted_josephson_ring_frequency_1(self):
        C = 1e-13
        L = 1e-8
        w,k,A,chi = self.shunted_josephson_ring_parameters(C,L)
        self.assertRelativelyClose(w[1],1/np.sqrt(L*C)/2./np.pi)

    def test_shunted_josephson_ring_anharmonicity_0(self):
        C = 1e-13
        L = 1e-8
        w,k,A,chi = self.shunted_josephson_ring_parameters(C,L)
        self.assertRelativelyClose(A[0],e**2/2./(8*C)/h)

    def test_shunted_josephson_ring_anharmonicity_1(self):
        C = 1e-13
        L = 1e-8
        w,k,A,chi = self.shunted_josephson_ring_parameters(C,L)
        self.assertRelativelyClose(A[1],e**2/2./(8*C)/h)

    '''
    Coupled transmon/RLC
    '''
    def test_sweeping_LJ_in_fkAchi(self):
        cir = core.Network([
            core.C(0,1,100e-15),
            core.J(0,1,'L_J'),
            core.C(1,2,1e-15),
            core.C(2,0,100e-15),
            core.L(2,0,10e-9),
            core.R(2,0,1e6)
            ])
        [cir.f_k_A_chi(L_J=x) for x in [1e-9,2e-9]]

class TestGraphics(TestCaseAppended):   

    def test_error_when_trying_to_plot_from_Network_show(self):
        circuit = core.Network([
            core.C(0,1,'C'),
            core.J(0,1,'Lj')
        ])
        with self.assertRaises(TypeError):
            circuit.show()
            
    def test_error_when_trying_to_plot_from_Network_show_normal_modes(self):
        circuit = core.Network([
            core.C(0,1,'C'),
            core.J(0,1,'Lj')
        ])
        with self.assertRaises(TypeError):
            circuit.show_normal_mode()
        

    def test_generate_graphics(self):
        import qucat.src._generate_graphics

class TestNetworkAnalysis(TestCaseAppended):

    def test_transfer_left_right_port_identical(self):
        '''
        Trivial cases
        '''
        net = core._Network([core.R(0,1,'Z2')])
        self.assertEqual(net.transfer(0,1,0,1),1)

    def test_transfer_left_right_port_indentical_inverted(self):
        '''
        Trivial cases
        '''
        net = core._Network([core.R(0,1,'Z2')])
        self.assertEqual(net.transfer(0,1,1,0),-1)

    def test_transfer_voltage_divider(self):
        '''
        Voltage divider, see:
        https://en.wikipedia.org/wiki/Voltage_divider
        We add an extra resistor between Vin and ground
        '''

        # Compute the bridge transfer function
        net = core._Network([
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
        net = core._Network([
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
        net = core._Network([
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
        net = core._Network([
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

    def test_open_or_series_check(self):
        with self.assertRaises(ValueError):
            net = core._Network([
                core.R(0,1,'Z'),
                core.C(1,2,'Z'),
                core.J(2,3,'Z'),
                ])

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