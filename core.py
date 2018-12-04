from __future__ import print_function
import lcapy
import sympy as sp
import numpy as np
from scipy.constants import e,pi,h
from sympy.core.mul import Mul,Pow,Add
from copy import deepcopy
import matplotlib.pyplot as plt
import math

class L(lcapy.oneport.L):
    def __init__(self, label):
        super(L,self).__init__(label)
        self.label = label
class J(lcapy.oneport.L):
    def __init__(self, label):
        super(J,self).__init__(label)
        self.label = label
class R(lcapy.oneport.R):
    def __init__(self, label):
        super(R, self).__init__(label)
        self.label = label
class C(lcapy.oneport.C):
    def __init__(self, label):
        super(C, self).__init__(label)
        self.label = label


def admittance(circuit):
    '''
    Given an Lcapy circuit, returns the admittance.
    '''
    if type(circuit) == lcapy.oneport.Par:
        return Add(admittance(circuit.args[0]),admittance(circuit.args[1]))
    elif type(circuit) == lcapy.oneport.Ser:
        return 1/Add(1/admittance(circuit.args[0]),1/admittance(circuit.args[1]))
    elif type(circuit) == L or type(circuit) == J:
        return -sp.I*Mul(1/sp.Symbol('w'),1/sp.Symbol(circuit.label,real=True,nonnegative = True))
    elif type(circuit) == C:
        return sp.I*Mul(sp.Symbol('w'),sp.Symbol(circuit.label,real=True,nonnegative = True))
    elif type(circuit) == R:
        return 1/sp.Symbol(circuit.label,real=True,nonnegative = True)

def remove_resistances(circuit):
    '''
    Given an Lcapy circuit, returns the same circuit without resistances
    '''
    if type(circuit) in [L,C,J]:
        return circuit
    elif type(circuit.args[0])==R:
        return circuit.args[1]
    elif type(circuit.args[1])==R:
        return circuit.args[0]
    elif type(circuit) == lcapy.oneport.Par:
        return lcapy.Par(remove_resistances(circuit.args[0]),remove_resistances(circuit.args[1]))
    elif type(circuit) == lcapy.oneport.Ser:
        return lcapy.Ser(remove_resistances(circuit.args[0]),remove_resistances(circuit.args[1]))

def get_all_circuit_elements(circuit):
    
    circuit_elements = {}
    
    def recursive_function(circuit):
        if type(circuit) == lcapy.oneport.Par:
            recursive_function(circuit.args[0])
            recursive_function(circuit.args[1])
        elif type(circuit) == lcapy.oneport.Ser:
            recursive_function(circuit.args[0])
            recursive_function(circuit.args[1])
        else:
            circuit_elements[circuit.label] = circuit
    
    recursive_function(circuit)
    return circuit_elements

def get_unreduced_node_circuit(circuit):
    # Classify the connections between nodes into 
        # wires (to be removed), capacitances and inductances
        wires = []
        capacitances = []
        inductances = []
        resistances = []

        for line in circuit.netlist().split('\n'):
            first_split = line.split((' '))
            
            connection = []
            if first_split[0]=='W':
                connection.append(int(first_split[1]))
                connection.append(int(first_split[2][:-1]))
                wires.append(connection)
            else:
                connection.append(int(first_split[1]))
                connection.append(int(first_split[2]))
                
                if '{' in line and '}' in line:
                    start = line.index('{') + 1
                    end = line.index( '}', start )
                    connection.append( line[start:end])
                else:
                    connection.append(first_split[3][:-1])
                
                if first_split[0]=='C':
                    capacitances.append(connection)
                elif first_split[0]=='L':
                    inductances.append(connection)
                elif first_split[0]=='R':
                    resistances.append(connection)
        return capacitances,inductances,resistances,wires

def get_node_circuit(circuit):
        capacitances,inductances,resistances,wires = get_unreduced_node_circuit(circuit)

        # Group nodes connected by wires into "chains" of nodes

        def merge_chains(chains,i,j):
            to_add = chains[j]
            del chains[j]
            chains[i] = chains[i]+to_add
            return chains


        chains = []
        for el in wires:
            added = False
            for i,ch in enumerate(chains):
                if (el[0] in ch) and (el[1] in ch):
                    added = True
                elif (el[0] in ch):
                    for j,ch2 in enumerate(chains):
                        if el[1] in ch2:
                            chains = merge_chains(chains,i,j)
                            added = True
                    if added == False:
                        ch.append(el[1])
                        added = True
                elif (el[1] in ch):
                    for j,ch2 in enumerate(chains):
                        if el[0] in ch2:
                            chains = merge_chains(chains,i,j)
                            added = True
                    if added == False:
                        ch.append(el[0])
                        added = True
            if added == False:
                chains.append(el)

        # Rewrite the capacitance and inductance 
        # conections in terms of a minimal number of nodes

        def reduce_node(node):
            for i,ch in enumerate(chains):
                if node in ch:
                    return i    

        for cpi in [capacitances,inductances,resistances]:
            for el in cpi:
                el[0] = reduce_node(el[0])
                el[1] = reduce_node(el[1])


        N_nodes = len(chains)
        return capacitances,inductances,resistances,N_nodes


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
        self.all_circuit_elements = get_all_circuit_elements(circuit)

        Y_numer = sp.numer(sp.together(Y))
        Y_poly = sp.collect(sp.expand(Y_numer),sp.Symbol('w'))
        self.Y_poly = Y_poly
        self.Y_poly_order = sp.polys.polytools.degree(Y_poly,gen = sp.Symbol('w'))
        Y_poly_coeffs = [Y_poly.coeff(sp.Symbol('w'),n) for n in range(self.Y_poly_order+1)[::-1]]
        self.Y_poly_num = sp.utilities.lambdify([elt for elt in self.all_circuit_elements.keys()],Y_poly_coeffs,"numpy")
        self.N_modes = int(self.Y_poly_order/2)
        self.dY = sp.diff(Y,sp.Symbol('w'))
        self.dY_num = sp.utilities.lambdify([elt for elt in self.all_circuit_elements.keys()+["w"]],self.dY,"numpy")
    
    def draw(self,
        y_total = 0.6,
        x_total=1.,
        y_C = 0.3,
        x_C = 0.2,
        x_J = 0.2,
        x_LR = 0.7,
        y_LR = 0.3,
        fontsize = 15,
        text_position = [0.,-0.35],
        plot = True,
        full_output = False):
        '''
        Draws the circuit
        '''
        line_type = []
                
        def draw(circuit, pos = [0.,0.], which = 'tc'):
            if type(circuit) == lcapy.oneport.Par:
                return draw_par(circuit.args[0],circuit.args[1],pos,which)
            elif type(circuit) == lcapy.oneport.Ser:
                return draw_ser(circuit.args[0],circuit.args[1],pos,which)
            elif type(circuit) == L:
                return draw_L(pos,which,circuit)
            elif type(circuit) == R:
                return draw_R(pos,which,circuit)
            elif type(circuit) == C:
                return draw_C(pos,which,circuit)
            elif type(circuit) == J:
                return draw_J(pos,which,circuit)


        def draw_par(el1,el2,pos,which):
            pos_x1,pos_y1,w1,h1,x1,y1,elt_names1 = draw(el1,[0.,0.],'t')
            pos_x2,pos_y2,w2,h2,x2,y2,elt_names2 = draw(el2,[0.,-h1],'t')

            w = max(w1,w2)
            h = h1+h2
            x = x1+x2
            y = y1+y2

            # Connect parallel elements

            x+=[np.array([-w/2.,-w/2.])]
            y+=[np.array([-h1/2.,-h1-h2/2.])]
            line_type.append('wire')

            x+=[np.array([+w/2.,+w/2.])]
            y+=[np.array([-h1/2.,-h1-h2/2.])]
            line_type.append('wire')

            if np.argmin([w1,w2]) == 0:
                x+=[np.array([-w/2.,-w1/2.])]
                y+=[np.array([-h1/2.,-h1/2.])]
                line_type.append('wire')
                x+=[np.array([+w/2.,+w1/2.])]
                y+=[np.array([-h1/2.,-h1/2.])]
                line_type.append('wire')
            if np.argmin([w1,w2]) == 1:
                x+=[np.array([-w/2.,-w2/2.])]
                y+=[np.array([-h1-h2/2.,-h1-h2/2.])]
                line_type.append('wire')
                x+=[np.array([+w/2.,+w2/2.])]
                y+=[np.array([-h1-h2/2.,-h1-h2/2.])]
                line_type.append('wire')

            # Move into position

            if which == 'l':
                return shift(pos_x1+pos_x2,pos[0]+w/2.),shift(pos_y1+pos_y2,pos[1]+h/2.),w,h,shift(x,pos[0]+w/2.),shift(y,pos[1]+h/2.),elt_names1+elt_names2
            else:
                return shift(pos_x1+pos_x2,pos[0]),shift(pos_y1+pos_y2,pos[1]),w,h,shift(x,pos[0]),shift(y,pos[1]),elt_names1+elt_names2


        def draw_ser(el1,el2,pos,which):
            pos_x1,pos_y1,w1,h1,x1,y1,elt_names1 = draw(el1,[0.,0.],'l')
            pos_x2,pos_y2,w2,h2,x2,y2,elt_names2 = draw(el2,[w1,0.],'l')
            w = w1+w2
            h = max(h1,h2)
            x = x1+x2
            y = y1+y2

            if which == 't':
                return shift(pos_x1+pos_x2,pos[0]-w/2.),shift(pos_y1+pos_y2,pos[1]-h/2.),w,h,shift(x,pos[0]-w/2.),shift(y,pos[1]-h/2.),elt_names1+elt_names2
            else:
                return shift(pos_x1+pos_x2,pos[0]),shift(pos_y1+pos_y2,pos[1]),w,h,shift(x,pos[0]),shift(y,pos[1]),elt_names1+elt_names2


        def shift(to_shift,shift):
            for i,_ in enumerate(to_shift):
                to_shift[i]+= shift
            return to_shift

        def draw_L(pos,which,circuit):
            x = np.linspace(0.5,6.,1001)
            y = -np.sin(2.*np.pi*x)
            x = np.cos(2.*np.pi*x)+2.*x
            line_type.append('L')
            return finish_drawing_LR([x],[y],pos,which,circuit)

        def draw_R(pos,which,circuit):
            x = np.linspace(-0.25,4.25,1001)
            height = 1.
            period = 1.
            a = height*2.*(-1.+2.*np.mod(np.floor(2.*x/period),2.))
            b = -height*2.*np.mod(np.floor(2.*x/period),2.)
            y = (2.*x/period - np.floor(2.*x/period))*a+b+height
            line_type.append('R')
            return finish_drawing_LR([x],[y],pos,which,circuit)


        def draw_C(pos,which,circuit):
            
            x = [
                np.array([0.,(x_total-x_C)/2.]),
                np.array([(x_total+x_C)/2.,x_total]),
                np.array([(x_total-x_C)/2.,(x_total-x_C)/2.]),
                np.array([(x_total+x_C)/2.,(x_total+x_C)/2.]),
            ]
            y = [
                np.array([0.,0.]),
                np.array([0.,0.]),
                np.array([-y_C/2.,y_C/2.]),
                np.array([-y_C/2.,y_C/2.]),
            ]
            line_type.append('wire')
            line_type.append('wire')
            line_type.append('C')
            line_type.append('C')
            
            for i,_ in enumerate(x):
                if which == 'l':
                    x[i]+=pos[0]
                    y[i]+=pos[1]
                    pos_x = pos[0]+x_total/2.
                    pos_y = pos[1]
                elif which == 't':
                    x[i]+=pos[0]-x_total/2.
                    y[i]+=pos[1]-y_total/2.
                    pos_x = pos[0]
                    pos_y = pos[1]-y_total/2.

            return [pos_x],[pos_y],x_total,y_total,x,y,[circuit.label]
        
        def draw_J(pos,which,circuit):
            
            x = [
                np.array([0.,(x_total-x_J)/2.]),
                np.array([(x_total+x_J)/2.,x_total]),
                np.array([(x_total-x_J)/2.,(x_total+x_J)/2.]),
                np.array([(x_total-x_J)/2.,(x_total+x_J)/2.]),
                np.array([(x_total-x_J)/2.,(x_total+x_J)/2.])
            ]
            y = [
                np.array([0.,0.]),
                np.array([0.,0.]),
                np.array([0.,0.]),
                np.array([-1.,1.])*x_J/2.,
                np.array([1.,-1.])*x_J/2.
            ]
            line_type.append('wire')
            line_type.append('wire')
            line_type.append('wire')
            line_type.append('J')
            line_type.append('J')
            
            for i,_ in enumerate(x):
                if which == 'l':
                    x[i]+=pos[0]
                    y[i]+=pos[1]
                    pos_x = pos[0]+x_total/2.
                    pos_y = pos[1]
                elif which == 't':
                    x[i]+=pos[0]-x_total/2.
                    y[i]+=pos[1]-y_total/2.
                    pos_x = pos[0]
                    pos_y = pos[1]-y_total/2.

            return [pos_x],[pos_y],x_total,y_total,x,y,[circuit.label]

        def finish_drawing_LR(x_list,y_list,pos,which,circuit):
            
            
            x_min = min([np.amin(xx) for xx in x_list])
            for i,_ in enumerate(x_list):
                x_list[i] -= x_min

            x_max = max([np.amax(xx) for xx in x_list])
            for i,_ in enumerate(x_list):
                x_list[i] *= x_LR/x_max
                x_list[i] +=(x_total-x_LR)/2.   

            x_min = min([np.amin(xx) for xx in x_list])
            x_max = max([np.amax(xx) for xx in x_list])
            x_list += [np.array([0.,x_min])]
            x_list += [np.array([x_max,x_total])]
            line_type.append('wire')
            line_type.append('wire')

            for i,x in enumerate(x_list):
                if which == 'l':
                    x_list[i]+=pos[0]
                elif which == 't':
                    x_list[i]+=pos[0]-x_total/2.


            y_max = max([np.amax(yy) for yy in y_list])
            for i,_ in enumerate(y_list):
                y_list[i] *=y_LR/2./y_max

            y_list += [np.array([0.,0.])]
            y_list += [np.array([0.,0.])]

            for i,y in enumerate(y_list):
                if which == 'l':
                    y_list[i]+=pos[1]
                elif which == 't':
                    y_list[i]+=pos[1]-y_total/2.

            if which == 'l':
                pos_x = pos[0]+x_total/2.
                pos_y = pos[1]
            elif which == 't':
                pos_x = pos[0]
                pos_y = pos[1]-y_total/2.

            return [pos_x],[pos_y],x_total,y_total,x_list,y_list,[circuit.label]
        
        element_x,element_y,w,h,xs,ys,element_names = draw(self.circuit,which='l')
        figsize_scaling = 1.5
        fig,ax = plt.subplots(1,1,figsize = (w*figsize_scaling,h*figsize_scaling))
        for i in range(len(xs)):
            if line_type[i] == 'wire':
                ax.plot(xs[i],ys[i],color = 'black',lw=1)
            elif line_type[i] == 'L':
                ax.plot(xs[i],ys[i],color = 'black',lw=2)
            elif line_type[i] == 'C':
                ax.plot(xs[i],ys[i],color = 'black',lw=8)
            elif line_type[i] == 'J':
                ax.plot(xs[i],ys[i],color = 'black',lw=8)
            elif line_type[i] == 'R':
                ax.plot(xs[i],ys[i],color = 'black',lw=2)
        for i,el in enumerate(element_names):
            plt.text(element_x[i]+text_position[0],element_y[i]+text_position[1],'$'+el+'$',fontsize=fontsize,ha='center')
        ax.set_axis_off()
        if plot:
            plt.tight_layout()
            plt.show()
        if full_output:
            return element_x,element_y,fig,ax

        
    def analytical_solution(self,simplify = False):
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

        facts = [sp.Q.positive(sp.Symbol(x)) for x in self.all_circuit_elements.keys()]
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
                    if simplify:
                        ws.append(sp.simplify(w))
                        As.append(sp.simplify(2*sp.Symbol('e')**2/sp.Symbol('h')*Mul(1/sp.Symbol(self.L_J),\
                                Mul(Pow(1/ImdY_lossless.subs({sp.Symbol('w'):w}),2),\
                                Pow(1/w,2)))))
                    else:
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
        args = [circuit_parameters[key] for key in self.all_circuit_elements.keys()]
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

    def get_L_and_C_matrix(self):

        '''
        Currently treats junctions as inductors
        '''

        cir_no_resistance = remove_resistances(self.circuit).simplify()
        capacitances,inductances,_,N_nodes = get_node_circuit(cir_no_resistance)

        # Write the upper part of the capacitance and 
        # inverse of the inductance matrix


        C_matrix = [[0 for col in range(N_nodes)] for row in range(N_nodes)]
        Lmin_matrix = [[0 for col in range(N_nodes)] for row in range(N_nodes)]

        for el in capacitances:
            nodes = [el[0],el[1]]
            C_matrix[min(nodes)][max(nodes)] = '-('+el[2]+')'

        for el in inductances:
            nodes = [el[0],el[1]]
            Lmin_matrix[min(nodes)][max(nodes)] = '-1/('+el[2]+')'

        print('''
            The diagonal elements are the opposite of sums of the row or column they belong to
            (excluding the lower triangle of the matrices).
            The matrices are symmetric so the upper triangle is identical to the lower triangle.
            We hence only give the upper triangle.
            Quantization can be performed by choosing a ground node and eliminating the row
            and column corresponding to that given node.
            ''')
        return C_matrix,Lmin_matrix


if __name__ == '__main__':
    
    ###LC circuit
    circuit = ((L('L1')|(L('L')+C('C')+J('L_J')+R('R')))+L('L'))|C('C')|J('L_J')+C('C')
    b = Bbox(circuit)
    print (b.get_L_and_C_matrix())