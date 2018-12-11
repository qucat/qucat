import sympy as sp
import numpy as np
from scipy.constants import e,pi,h,hbar
from sympy.core.mul import Mul,Pow,Add
from copy import deepcopy
id2 = sp.Matrix([[1,0],[0,1]])

class Circuit(object):
    """docstring for Circuit"""
    def __init__(self):
        self.id_ = None 
        self.parent = None
        self.head = None 

    def is_equal_to(self,other_circuit):
        return self.id_ == other_circuit.id_

    def __add__(self, other_circuit):
        return Series(self, other_circuit)

    def __or__(self, other_circuit):
        return Parallel(self, other_circuit)


class Connection(Circuit):
    """docstring for Connection"""
    def __init__(self, left, right):
        super(Connection, self).__init__()
        self.left = left
        self.right = right
        self.set_parenthood()

    def set_parenthood(self):
        self.left.parent = self
        self.right.parent = self

    def remove_resistances(self):
        if type(self.left)==R:
            return self.right
        elif type(circuit.right)==R:
            return self.left
        else:
            return Connection(self.left,self.right)


    def set_headConnection_componentIds_componentDict(self,head = None):
        if head == None:
            self.head = self
            self.id_counter = 0
            self.component_dict = {}
            self.inductors = []
            self.capacitors = []
            self.junctions = []
            self.resistors = []
        else:
            self.head = head

        self.id_ = self.head.id_counter
        self.head.id_counter += 1
        self.left.set_headConnection_componentIds_componentDict(self.head)
        self.right.set_headConnection_componentIds_componentDict(self.head)

class Series(Connection):
    """docstring for Series"""
    def __init__(self, left,right):
        super(Series, self).__init__(left,right)

    def admittance(self):
        return 1/Add(1/self.left.admittance(),
            1/self.right.admittance())

    def set_flux_transforms(self,ABCD = id2):
        self.right.set_flux_transforms(
            ABCD*sp.Matrix([[1,1/self.left.admittance()],[0,1]]))
        self.left.set_flux_transforms(
            ABCD*sp.Matrix([[1,1/self.right.admittance()],[0,1]]))

class Parallel(Connection):
    """docstring for Parallel"""
    def __init__(self, left,right):
        super(Parallel, self).__init__(left,right)

        #########################
        # Head connection variables
        #########################
        self.Q_min = 1.
        self.circuit_rotated = None
        self.Y = None
        self.Y_poly_coeffs = None
        self.w_cpx = None
        self.flux_transforms_set = False
        self.dY = None

        
    def admittance(self):
        return Add(self.left.admittance(),
            self.right.admittance())

    def set_flux_transforms(self,ABCD = id2):
        self.right.set_flux_transforms(
            ABCD*sp.Matrix([[1,0],[self.left.admittance(),1]]))
        self.left.set_flux_transforms(
            ABCD*sp.Matrix([[1,0],[self.right.admittance(),1]]))



    #########################
    # Head connection functions
    #########################

    def set_circuit_rotated(self):
        if self.circuit_rotated == None:

            if self.id_ == None:
                self.set_headConnection_componentIds_componentDict()

            if len(self.junctions)>0:
                self.circuit_rotated = rotate(self.junctions[0]) 
            else:
                self.circuit_rotated = rotate_circuit(self.circuit,self.inductors[0]) 


    def set_Y(self):
        if self.Y == None:
            self.set_circuit_rotated()
            self.Y = self.circuit_rotated.admittance() 
    
    def set_Y_poly_coeffs(self):
        if self.Y_poly_coeffs == None:
            self.set_Y()
            Y_numer = sp.numer(sp.together(self.Y) )       # Extract the numerator of Y(w)
            Y_poly = sp.collect(sp.expand(Y_numer),sp.Symbol('w')) # Write numerator as polynomial in omega
            Y_poly_order = sp.polys.polytools.degree(Y_poly,gen = sp.Symbol('w')) # Order of the polynomial
            self.Y_poly_coeffs = [complex(Y_poly.coeff(sp.Symbol('w'),n)) for n in range(Y_poly_order+1)[::-1]] # Get polynomial coefficients

    
    def set_w_cpx(self):
        if self.w_cpx == None:
            self.set_Y_poly_coeffs()
            ws_cpx = np.roots(self.Y_poly_coeffs)

            # take only roots with a positive real part (i.e. freq) 
            # and significant Q factors
            relevant_sols = np.argwhere((np.real(ws_cpx)>=0.)&(np.real(ws_cpx)>self.Q_min*np.imag(ws_cpx)))
            ws_cpx=ws_cpx[relevant_sols][:,0]
        
            # Sort solutions with increasing frequency
            self.w_cpx = np.sort(np.real(ws_cpx))
    
    def set_dY(self):
        if self.dY == None:
            self.set_Y()
            self.dY = sp.utilities.lambdify('w',sp.diff(self.Y,sp.Symbol('w')),'numpy')

    def compute_all_flux_transformations(self):
        if self.flux_transforms_set == False:
            self.circuit_rotated.set_flux_transforms()
            self.flux_transforms_set = True

    def eigenfrequencies(self):
        self.set_w_cpx()
        return np.real(self.w_cpx)/2./pi

    def loss_rates(self):
        return np.imag(self.w_cpx)/2./pi

    def Qs(self):
        return np.real(self.w_cpx)/np.imag(self.w_cpx)

    def anharmonicities(self):
        self.set_dY()
        if len(self.junctions) == 0:
            print "There are no junctions and hence no anharmonicity in the circuit"

        elif len(self.junctions) == 1:
            self.junctions[0].flux_wr_ref = lambda w: 1.
            return np.absolute(self.junctions[0].anharmonicity(self.w_cpx))/h

        else:
            self.compute_all_flux_transformations()
            return sum([np.absolute(j.anharmonicity(self.w_cpx)) for j in self.junctions])/h


class Component(Circuit):
    """docstring for Component"""
    def __init__(self, value):
        super(Component, self).__init__()
        self.value = value
        self.flux = lambda x: 0.

    def remove_resistances(self):
        return self

    def set_headConnection_componentIds_componentDict(self,head = None):

        self.head = head
        self.id_ = self.head.id_counter
        self.head.id_counter += 1
        self.head.component_dict[self.id_] = self
        self.add_component_category_list()

    def set_flux_transforms(self,ABCD = id2):
        ABCD = ABCD*sp.Matrix([[1,0],[self.admittance(),1]])
        self.flux_wr_ref = sp.utilities.lambdify(['w'],1/ABCD[0,0],"numpy")

    def flux(self,w):
        ImdY = np.imag(self.head.dY(w))
        return self.flux_wr_ref(w)*sp.sqrt(hbar/w/ImdY)

    def charge(self,w):
        return hbar/2./self.flux(w)

    def voltage(self,w):
        return self.flux(w)*1j*w

    def current(self,w):
        Y = sp.utilities.lambdify(['w'],self.admittance,"numpy")
        return self.voltage(w)*Y(w)

class L(Component):
    def __init__(self, value):
        super(L,self).__init__(value)
    def admittance(self):
        return -sp.I*Mul(1/sp.Symbol('w'),1/self.value)

    def add_component_category_list(self):
        self.head.inductors.append(self)

class J(L):
    def __init__(self, value):
        super(L,self).__init__(value)

    def add_component_category_list(self):
        self.head.junctions.append(self)

    def anharmonicity(self,w):
        ImdY = np.imag(self.head.dY(w))
        return self.flux_wr_ref(w)**4*2.*e**2/self.value/w**2/ImdY**2
        
class R(Component):
    def __init__(self, value):
        super(R, self).__init__(value)
    def admittance(self):
        return 1/self.value

    def add_component_category_list(self):
        self.head.resistors.append(self)

class C(Component):
    def __init__(self, value):
        super(C, self).__init__(value)
    def admittance(self):
        return sp.I*Mul(sp.Symbol('w'),self.value)

    def add_component_category_list(self):
        self.head.capacitors.append(self)


def rotate(circuit_element):
    
    # Check if the tree is already correctly rotated
    if circuit_element.parent.parent == None:
        return circuit_element.parent
    
    def recursive_function(elt):
        parent_connection = elt.parent
        if parent_connection.parent == None:
            if parent_connection.left.is_equal_to(elt):
                return parent_connection.right
            elif parent_connection.right.is_equal_to(elt):
                return parent_connection.left
        else:
            if parent_connection.left.is_equal_to(elt):
                parent_connection.left = recursive_function(parent_connection)
            elif parent_connection.right.is_equal_to(elt):
                parent_connection.right = recursive_function(parent_connection)
            return parent_connection
        
    rotated_circuit = circuit_element|recursive_function(circuit_element)
    rotated_circuit.id_ = 0
    
    def reset_all_parenthoods(elt):
        if type(elt) == Connection:
            elt.set_parenthood()
            set_parenthoods(elt.left)
            set_parenthoods(elt.right)
    
    reset_all_parenthoods(rotated_circuit)

    return rotated_circuit

if __name__ == '__main__':
    for i in range(1):
        circuit = (J(10e-9)|((C(100e-15))|(C(1e-15)+(C(100e-15)|L(10e-9)|R(1e8)))))|(J(10e-9)|((C(100e-15))|(C(1e-15)+(C(50e-15)|L(10e-9)|R(1e8)))))+(J(10e-9)|((C(100e-15))|(C(1e-15)+(C(100e-15)|L(10e-9)|R(1e8)))))
        circuit.eigenfrequencies()
        circuit.anharmonicities()