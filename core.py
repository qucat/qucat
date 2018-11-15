from __future__ import print_function
from lcapy import L,C,R
import lcapy
import sympy as sp
import numpy as np
from scipy.constants import e,pi,h
from sympy.core.mul import Mul,Pow,Add
from copy import deepcopy
import matplotlib.pyplot as plt
import math

def admittance(circuit):
    '''
    Given an Lcapy circuit, returns the admittance.
    '''
    if type(circuit) == lcapy.oneport.Par:
        return Add(admittance(circuit.args[0]),admittance(circuit.args[1]))
    elif type(circuit) == lcapy.oneport.Ser:
        return 1/Add(1/admittance(circuit.args[0]),1/admittance(circuit.args[1]))
    elif type(circuit) == lcapy.oneport.L:
        return -sp.I*Mul(1/sp.Symbol('w'),1/sp.Symbol(circuit.args[0],real=True))
    elif type(circuit) == lcapy.oneport.C:
        return sp.I*Mul(sp.Symbol('w'),sp.Symbol(circuit.args[0],real=True))
    elif type(circuit) == lcapy.oneport.R:
        return 1/sp.Symbol(circuit.args[0],real=True)

def remove_resistances(circuit):
    '''
    Given an Lcapy circuit, returns the same circuit without resistances
    '''
    if type(circuit) == lcapy.oneport.L:
        return circuit
    elif type(circuit) == lcapy.oneport.C:
        return circuit
    elif type(circuit.args[0])==lcapy.oneport.R:
        return circuit.args[1]
    elif type(circuit.args[1])==lcapy.oneport.R:
        return circuit.args[0]
    elif type(circuit) == lcapy.oneport.Par:
        return lcapy.Par(remove_resistances(circuit.args[0]),remove_resistances(circuit.args[1]))
    elif type(circuit) == lcapy.oneport.Ser:
        return lcapy.Ser(remove_resistances(circuit.args[0]),remove_resistances(circuit.args[1]))

class Bbox(object):
    '''
    Given a circuit, the Bbox allows one to draw the circuit, and compute numerically
    and analytically (not yet) the frequency, dissipation rate and anharmonicity
    of the circuits different modes.
    '''
    def __init__(self, circuit, L_J = "L_J", Q_min = 1.):
        self.circuit = circuit#.simplify()
        self.L_J = L_J
        self.Q_min = Q_min
        
        Y = admittance(self.circuit)
        self.all_circuit_elements = [str(x) for x in list(Y.free_symbols)]
        self.all_circuit_elements.remove('w')

        Y_numer = sp.numer(sp.together(Y))
        Y_poly = sp.collect(sp.expand(Y_numer),sp.Symbol('w'))
        self.Y_poly = Y_poly
        self.Y_poly_order = sp.polys.polytools.degree(Y_poly,gen = sp.Symbol('w'))
        Y_poly_coeffs = [Y_poly.coeff(sp.Symbol('w'),n) for n in range(self.Y_poly_order+1)[::-1]]
        self.Y_poly_num = sp.utilities.lambdify([elt for elt in self.all_circuit_elements],Y_poly_coeffs,"numpy")
        self.N_modes = int(self.Y_poly_order/2)
        self.dY = sp.diff(Y,sp.Symbol('w'))
        self.dY_num = sp.utilities.lambdify([elt for elt in self.all_circuit_elements+["w"]],self.dY,"numpy")
    
    def draw(self):
        '''
        Draws the circuit
        '''
        self.circuit.draw()
        plt.show()
        
    def analytical_solution(self):
        '''
        Attempts to return analytical expressions for the frequency and anharmonicity
        Does not attempt to calculate the dissipation, since sympy does not simplify
        real/imaginary parts of complex expressions well
        '''
        self.circuit_lossless = remove_resistances(self.circuit)
        Y_lossless = admittance(self.circuit_lossless)
        Y_numer_lossless = sp.numer(sp.together(Y_lossless))
        Y_poly_lossless = sp.collect(sp.expand(Y_numer_lossless),sp.Symbol('w'))
        ImdY_lossless = sp.diff(Y_lossless,sp.Symbol('w')).subs({sp.I:1})

        facts = [sp.Q.positive(sp.Symbol(x)) for x in self.all_circuit_elements]
        with sp.assuming(*facts):
        
            # Try and calculate analytical eigenfrequencies
            w_analytical = sp.solve(Y_poly_lossless,sp.Symbol('w'))

            # Check the number of solutions
            if len(w_analytical)==0:
                print ("No analytical solutions")
                return None

            ws = []
            As = []
            for w in w_analytical:
                w_num = w.evalf(subs={i:1. for i in w.free_symbols})
                if w_num>0:
                    ws.append(w)
                    As.append(2*sp.Symbol('e')**2/sp.Symbol('h')*Mul(1/sp.Symbol(self.L_J),\
                            Mul(Pow(1/ImdY_lossless.subs({sp.Symbol('w'):w}),2),\
                            Pow(1/w,2))))
            return ws,As
    
    def fkA(self,circuit_parameters):
        '''
        Input: dictionnary of circuit parameters
        Returns resonance frequencies (in Hz)
                dissipation rates (Hz)
                Anharmonicities (Hz)
        '''
        args = [circuit_parameters[key] for key in self.all_circuit_elements]
        ws_cpx = np.roots(self.Y_poly_num(*args))

        # take only roots with a positive real part (i.e. freq) 
        # and significant Q factors
        positive_sols = np.argwhere((np.real(ws_cpx)>=0.)&(np.real(ws_cpx)>self.Q_min*np.imag(ws_cpx)))
        ws_cpx=ws_cpx[positive_sols][:,0]
        
        ws = np.real(ws_cpx)

        # Sort solutions with increasing frequency
        increasing_frequencies = np.argsort(ws)
        ws = ws[increasing_frequencies]
        ks = np.imag(ws_cpx)[increasing_frequencies]
        
        
        args.append(ws)
        ImdY = np.imag(self.dY_num(*args))
        As = 2.*e**2/h/circuit_parameters[self.L_J]/ws**2/ImdY**2
        return np.concatenate(([ws/2./pi],[ks/2./pi],[As]))
    
    def normalmodes(self,circuit_parameters):
        '''
        Input:  dictionnary of circuit parameters
                one of these parameters may be given as a list or numpy array
        Returns resonance frequencies (in Hz)
                dissipation rates (Hz)
                Anharmonicities (Hz)
        '''
        N_lists = 0
        for key in circuit_parameters:
            if type(circuit_parameters[key]) == np.ndarray or type(circuit_parameters[key])== list:
                N_lists +=1
                key_list = key
        
        if N_lists == 0:
            return self.fkA(circuit_parameters)
        elif N_lists == 1:
            to_iterate = deepcopy(circuit_parameters[key_list])
            circuit_parameters[key_list] = to_iterate[0]
            ret = np.array([self.fkA(circuit_parameters)])
            for key_value in to_iterate[1:]:
                circuit_parameters[key_list] = key_value
                ret = np.concatenate((ret,[self.fkA(circuit_parameters)]),axis = 0)
            return np.transpose(ret,(2,1,0))
        else:
            "Can only iterate on one variable"

if __name__ == '__main__':
    
    # LC circuit
    
    b = Bbox(L('L_J') | C('C'))
    print(b.analytical_solution())


    # Typical cQED setup

    # b_cQED = Bbox(L('L_J') | C('C')| R('R_J') | (C('Cc')+(C('Cr')|L('Lr')|(C('Cf')+R('R_50')))))
    # flux = np.linspace(0.,0.4,103)
    # L_J_list = 7e-9/abs(np.cos(np.pi*flux))
    # to_plot =  np.array(b_cQED.normalmodes({
    # 'L_J':L_J_list,
    # 'C':100e-15,
    # 'Cc':10e-15,
    # 'Cf':1e-15,
    # 'R_50':50.,
    # 'Lr':10e-9,
    # 'Cr':100e-15,
    # 'R_J':1e6,
    # 'R_r':1e5
    # }))
    # fig,axarr = plt.subplots(3,1,figsize=(6,9),sharex=True)
    # for i,ax in enumerate(axarr):
    #     for j in range(b_cQED.N_modes):
    #         axarr[i].plot(flux,to_plot[j,i])
    # axarr[0].set_ylabel("Mode frequency")
    # axarr[1].set_ylabel("Mode dissipation")
    # axarr[2].set_ylabel("Mode anharmonicity")
    # axarr[2].set_xlabel("$\phi/\phi_0$")
    # plt.show()