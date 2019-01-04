import sympy as sp
import numpy as np
from scipy.constants import e,pi,h,hbar
from sympy.core.mul import Mul,Pow,Add
from copy import deepcopy
from json import load
import matplotlib.pyplot as plt
from numbers import Number
from math import floor
phi_0 = hbar/2./e
id2 = sp.Matrix([[1,0],[0,1]])
exponent_to_letter = {
    -18:'a',
    -15:'f',
    -12:'p',
    -9:'n',
    -6:r'$\mu$',
    -3:'m',
    0:'',
    3:'k',
    6:'M',
    9:'G',
    12:'T'
}

class BBQcircuit(object):
    """docstring for BBQcircuit"""
    def __init__(self, elements):
        super(BBQcircuit, self).__init__()
        self.elements = elements
        self.network = Network(elements)


        self.inductors = []
        self.capacitors = []
        self.junctions = []
        self.resistors = []
        self.no_value_components = []
        for elt in elements:
            elt.head = self
            elt.set_component_lists()

        if len(self.junctions)>1:
            self.ref_elt = self.junctions[0]
        else:
            self.ref_elt = self.inductors[0]

        self.Q_min = 1.
        self.Y = self.network.admittance(self.ref_elt.node_minus,self.ref_elt.node_plus) 
        Y_together = sp.together(self.Y)
        Y_numer = sp.numer(Y_together)       # Extract the numerator of Y(w)
        Y_poly = sp.collect(sp.expand(Y_numer),sp.Symbol('w')) # Write numerator as polynomial in omega
        Y_poly_order = sp.polys.polytools.degree(Y_poly,gen = sp.Symbol('w')) # Order of the polynomial
        self.Y_poly_coeffs_analytical = [Y_poly.coeff(sp.Symbol('w'),n) for n in range(Y_poly_order+1)[::-1]] # Get polynomial coefficients
        self.Y_poly_coeffs_analytical = [sp.utilities.lambdify(self.no_value_components,c,'numpy') for c in self.Y_poly_coeffs_analytical]

        self.dY = sp.utilities.lambdify(['w']+self.no_value_components,sp.diff(self.Y,sp.Symbol('w')),'numpy')

        for elt in elements:
            tr = self.network.transfer(self.ref_elt.node_minus,self.ref_elt.node_plus,elt.node_minus,elt.node_plus)
            elt.flux_wr_ref = sp.utilities.lambdify(['w']+self.head.no_value_components,tr,"numpy")

    def Y_poly_coeffs(self,**kwargs):
        return [complex(coeff(**kwargs)) for coeff in self.Y_poly_coeffs_analytical]

    def w_k_A_chi(self,pretty_print = False,**kwargs):
        
        list_element = None
        list_values = None
        for el,value in kwargs.items():
            try:
                iter(value)
            except TypeError:
                iterable = False
            else:
                iterable = True

            if iterable and list_element is None:
                list_element = el 
                list_values = value
            elif iterable and list_element is not None:
                raise ValueError("You can only iterate over the value of one element.")

        if pretty_print == True and list_element is not None:
            raise ValueError("Cannot pretty print since $%s$ does not have a unique value"%list_element)

        if list_element is None:

            to_return =  self.eigenfrequencies(**kwargs),\
                    self.loss_rates(**kwargs),\
                    self.anharmonicities(**kwargs),\
                    self.kerr(**kwargs)

            if pretty_print:
                N_modes = len(to_return[0])
                table_line = ""
                for i in range(4):
                    table_line += " %7s |"
                table_line += "\n"

                to_print = table_line%("mode"," freq. "," diss. "," anha. ")
                for i,w in enumerate(to_return[0]):
                    to_print+=table_line%tuple([str(i)]+[pretty_value(to_return[j][i])+'Hz' for j in range(3)])

                to_print += "\nKerr coefficients\n(diagonal = Kerr, off-diagonal = cross-Kerr)\n"


                table_line = ""
                for i in range(N_modes+1):
                    table_line += " %7s |"
                table_line += "\n"

                to_print += table_line%tuple(['mode']+[str(i)+'   ' for i in range(N_modes)])

                for i in range(N_modes):
                    line_elements = [str(i)]
                    for j in range(N_modes):
                        if i>=j:
                            line_elements.append(pretty_value(to_return[3][i][j])+'Hz')
                        else:
                            line_elements.append("")
                    to_print += table_line%tuple(line_elements)
                print(to_print)

            return to_return

        else:
            w = []
            k = []
            A = []
            kerr = []
            for value in list_values:
                kwargs_single = deepcopy(kwargs)
                kwargs_single[list_element] = value
                w.append(self.eigenfrequencies(**kwargs_single))
                k.append(self.loss_rates(**kwargs_single))
                A.append(self.anharmonicities(**kwargs_single))
                kerr.append(self.kerr(**kwargs_single))
            w = np.moveaxis(np.array(w),0,-1)
            k = np.moveaxis(np.array(k),0,-1)
            A = np.moveaxis(np.array(A),0,-1)
            kerr = np.moveaxis(np.array(kerr),0,-1)
            return w,k,A,kerr

    def check_kwargs(self,**kwargs):
        for key in kwargs:
            if key in self.no_value_components:
                pass
            elif key in [c.label for _,c in self.component_dict.iteritems()]:
                raise ValueError('The value of %s was already specified when constructing the circuit'%key)
            else:
                raise ValueError('%s is not the label of a circuit element'%key)

        for label in self.no_value_components:
            try:
                kwargs[label]
            except Exception as e:
                raise ValueError('The value of %s should be specified with the keyword argument %s=... '%(label,label))

    def set_w_cpx(self,**kwargs):
        self.check_kwargs(**kwargs)
        ws_cpx = np.roots(self.Y_poly_coeffs(**kwargs))

        # take only roots with a positive real part (i.e. freq) 
        # and significant Q factors
        relevant_sols = np.argwhere((np.real(ws_cpx)>=0.)&(np.real(ws_cpx)>self.Q_min*np.imag(ws_cpx)))
        ws_cpx=ws_cpx[relevant_sols][:,0]
    
        # Sort solutions with increasing frequency
        order = np.argsort(np.real(ws_cpx))
        self.w_cpx = ws_cpx[order]
    
    def anharmonicities_per_junction(self,pretty_print =False,**kwargs):
        self.set_w_cpx(**kwargs)

        if len(self.junctions) == 0:
            raise UserWarning("There are no junctions and hence no anharmonicity in the circuit")
            return []

        elif len(self.junctions) == 1:
            def flux_wr_ref(w,**kwargs):
                return 1.
            self.junctions[0].flux_wr_ref = flux_wr_ref
            return [self.junctions[0].anharmonicity(self.w_cpx,**kwargs)/h]

        else:
            self.compute_all_flux_transformations()
            return [j.anharmonicity(self.w_cpx,**kwargs)/h for j in self.junctions]

    def kerr(self,**kwargs):
        As =  self.anharmonicities_per_junction(**kwargs)
        N_modes = len(As[0])
        N_junctions = len(As)

        Ks = np.zeros((N_modes,N_modes))
        for i in range(N_modes):
            line = []
            for j in range(N_modes):
                for k in range(N_junctions):
                    if i==j:
                        Ks[i,j]+=np.absolute(As[k][i])
                    else:
                        Ks[i,j]+=2.*np.sqrt(np.absolute(As[k][i])*np.absolute(As[k][j]))
        return Ks

    def anharmonicities(self,**kwargs):
        Ks = self.kerr(**kwargs)
        return [Ks[i,i] for i in range(Ks.shape[0])]


class Network(object):
    """docstring for Network"""
    def __init__(self, netlist):

        # construct node_names
        self.nodes = []
        for elt in netlist:
            nodes = [elt.node_minus,elt.node_plus]
            for node in nodes:
                if node not in self.nodes:
                    self.nodes.append(node)

        # initialize netlist_dict
        self.net_dict = {}
        for node in self.nodes:
            self.net_dict[node] = {}
        for elt in netlist:
            self.connect(elt,elt.node_minus,elt.node_plus)

    def connect(self,element,node_minus,node_plus):
        try:
            self.net_dict[node_minus][node_plus] = self.net_dict[node_minus][node_plus]|element
        except KeyError:
            self.net_dict[node_minus][node_plus] = element

        try:
            self.net_dict[node_plus][node_minus] = self.net_dict[node_plus][node_minus]|element
        except KeyError:
            self.net_dict[node_plus][node_minus] = element

    def remove_node(self,node):

        connections = [x for x in self.net_dict[node].items()]

        # Compute sum of admittances
        sum_Y = sum([elt.admittance() for _,elt in connections])

        # Add mesh
        for i,(node_A,elt_A) in enumerate(connections):
            for node_B,elt_B in connections[i:]:
                Y = sum_Y/elt_A.admittance()/elt_B.admittance()
                self.connect(Admittance(Y),node_A,node_B)

        # Remove star
        for other_node in self.net_dict[node]:
            del self.net_dict[other_node]
        del self.net_dict[node]

    def admittance(self,node_minus,node_plus):
        network_to_reduce = deepcopy(self)
        for node in self.nodes:
            if node not in [node_minus,node_plus]:
                network_to_reduce.remove_node(node)
        return network_to_reduce.net_dict[node_minus][node_plus].admittance()

    def transfer(self,node_minus_minus,node_minus_plus,node_plus_minus,node_plus_plus):

        # Reduce network
        network_to_reduce = deepcopy(self)
        for node in self.nodes:
            if node not in [node_minus_minus,node_minus_plus,node_plus_minus,node_plus_plus]:
                network_to_reduce.remove_node(node)

        # Compute ABCD of lattice network
        # see https://www.globalspec.com/reference/71734/203279/10-11-lattice-networks
        # Network Analysis & Circuit (By M. Arshad )section 10.11: LATTICE NETWORKS
        Za = 1/network_to_reduce.net_dict[node_minus_plus][node_plus_plus].admittance()
        Zb = 1/network_to_reduce.net_dict[node_minus_minus][node_plus_plus].admittance()
        Zc = 1/network_to_reduce.net_dict[node_minus_plus][node_plus_minus].admittance()
        Zd = 1/network_to_reduce.net_dict[node_minus_minus][node_plus_minus].admittance()
        sum_Z = sum([Za,Zb,Zc,Zd])
        Z11 = (Za+Zb)*(Zd+Zc)/sum_Z
        Z21 = (Zb*Zc-Za*Zd)/sum_Z
        Z22 = (Za+Zc)*(Zd+Zb)/sum_Z

        # see Pozar
        ABCD = sp.Matrix([[
            Z11/Z21,
            Z11*Z22/Z21-Z21],[
            1/Z21,
            Z22/Z21
            ]])

        # Connect missing two elements
        Y_L = network_to_reduce.net_dict[node_minus_plus][node_minus_minus].admittance()
        Y_R = network_to_reduce.net_dict[node_plus_plus][node_plus_minus].admittance()
        ABCD_L = sp.Matrix([[1,0],[Y_L,1]])
        ABCD_R = sp.Matrix([[1,0],[Y_R,1]])
        ABCD = ABCD_L*ABCD*ABCD_R

        return 1/ABCD[0,0]



class Circuit(object):
    """docstring for Circuit"""
    def __init__(self,node_minus,node_plus):
        self.node_minus = node_minus
        self.node_plus = node_plus
        self.head = None 

    def __or__(self, other_circuit):
        return Parallel(self, other_circuit)


class Parallel(Circuit):
    """docstring for Connection"""
    def __init__(self, left, right):
        super(Parallel, self).__init__(node_minus = None,node_plus = None)

        # sets the two children circuit elements
        self.left = left
        self.right = right
   
    def admittance(self):
        return Add(
            self.left.admittance(),
            self.right.admittance())


class Component(Circuit):
    """docstring for Component"""
    def __init__(self, node_minus, node_plus, arg1 = None, arg2 = None):
        super(Component, self).__init__(node_minus, node_plus)
        self.label = None
        self.value = None

        if arg1 is None and arg2 is None:
            raise ValueError("Specify either a value or a label")
        for a in [arg1,arg2]:
            if a is None:
                pass
            elif type(a) is str:
                self.label = a
            else:
                self.value = float(a)

    def get_value(self,**kwargs):
        if self.value is not None:
            return self.value
        elif self.value is None and kwargs is not None:
            if self.label in [k for k in kwargs]:
                return kwargs[self.label]
        
        return sp.Symbol(self.label)

    def set_component_lists(self):
        pass


    def flux(self,w,**kwargs):
        ImdY = np.imag(self.head.dY(w,**kwargs))
        return complex(self.flux_wr_ref(w,**kwargs)*sp.sqrt(hbar/w/ImdY))

    def voltage(self,w,**kwargs):
        return complex(self.flux(w,**kwargs)*1j*w)

    def current(self,w,**kwargs):
        kwargs['w']=w
        Y = self.admittance()
        if isinstance(Y, Number):
            pass
        else:
            Y = Y.evalf(subs=kwargs)
        return complex(self.voltage(**kwargs)*Y)

    def charge(self,w,**kwargs):
        return self.current(w,**kwargs)/w

    def to_string(self):

        if self.label is None:
            return pretty_value(self.value)+self.unit
        elif self.value is None:
            return ("$%s$"%(self.label))
        else:
            return ("$%s = $"%(self.label))+pretty_value(self.value)+self.unit

class L(Component):
    def __init__(self,node_minus, node_plus, arg1 = None, arg2 = None):
        super(L,self).__init__(node_minus, node_plus,arg1,arg2)
        self.unit = 'H'
    def admittance(self):
        return -sp.I*Mul(1/sp.Symbol('w'),1/self.get_value())

    def set_component_lists(self):
        super(L, self).set_component_lists()
        self.head.inductors.append(self)


class J(L):
    def __init__(self, node_minus, node_plus,arg1 = None, arg2 = None,use_E=False,use_I=False):
        super(J,self).__init__(node_minus, node_plus,arg1,arg2)

        self.use_E = use_E
        self.use_I = use_I
        if self.use_E:
            self.unit = 'Hz'
        elif self.use_I:
            self.unit = 'A'
        else:
            self.unit = 'H'

    def get_value(self,**kwargs):
        value = super(J,self).get_value(**kwargs)
        if (self.use_E == False) and (self.use_I ==False):
            return value
        elif (self.use_E == True) and (self.use_I ==False):
            L = (hbar/2./e)**2/(value*h) # E is assumed to be provided in Hz
            return L
        elif (use_E == False) and (use_I ==True):
            L = (hbar/2./e)/value
            return L
        else:
            raise ValueError("Cannot set both use_E and use_I to True")

    def set_component_lists(self):
        super(J, self).set_component_lists()
        self.head.junctions.append(self)

    def anharmonicity(self,w,**kwargs):
        ImdY = np.imag(self.head.dY(w,**kwargs))
        return self.flux_wr_ref(w,**kwargs)**4*2.*e**2/self.get_value(**kwargs)/w**2/ImdY**2


class R(Component):
    def __init__(self,node_minus, node_plus, arg1 = None, arg2 = None):
        super(R, self).__init__(node_minus, node_plus,arg1,arg2)
        self.unit = r'$\Omega$'

    def admittance(self):
        return 1/self.get_value()

    def set_component_lists(self):
        super(R, self).set_component_lists()
        self.head.resistors.append(self)


class C(Component):
    def __init__(self, node_minus, node_plus,arg1 = None, arg2 = None):
        super(C, self).__init__(node_minus, node_plus,arg1,arg2)
        self.unit = 'F'
    def admittance(self):
        return sp.I*Mul(sp.Symbol('w'),self.get_value())

    def set_component_lists(self):
        super(C, self).set_component_lists()
        self.head.capacitors.append(self)

class Admittance(Component):
    def __init__(self, Y):
        super(Admittance, self).__init__(node_minus = None, node_plus = None)
        self.Y = Y

    def admittance(self):
        return self.Y


def pretty_value(v,use_power_10 = False):
    exponent = floor(np.log10(v))
    exponent_3 = exponent-(exponent%3)
    float_part = v/(10**exponent_3)
    if use_power_10 or v>=1e15 or v<1e-18:
        if exponent_3 == 0:
            exponent_part = ''
        else:
            exponent_part = r'$\times 10^{%d}$'%exponent_3
    else:
        exponent_part = ' '+exponent_to_letter[exponent_3]
    if float_part>=10.:
        pretty = "%.0f%s"%(float_part,exponent_part)
    else:
        pretty = "%.1f%s"%(float_part,exponent_part)
    return pretty

def check_there_are_no_iterables_in_kwarg(**kwargs):
    for el,value in kwargs.items():
        try:
            iter(value)
        except TypeError:
            pass
        else:
            raise ValueError("This function accepts no lists or iterables as input.")

if __name__ == '__main__':
    cQED_circuit = BBQcircuit([
        C(0,1,100e-15),
        J(0,1,1e-9)])
    cQED_circuit