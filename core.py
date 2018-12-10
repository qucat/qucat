from __future__ import print_function
import lcapy
import sympy as sp
import numpy as np
from scipy.constants import e,pi,h,hbar
phi_0 = hbar/2./e
from sympy.core.mul import Mul,Pow,Add
from copy import deepcopy
import matplotlib.pyplot as plt
from rotations import rotate_circuit
from circuit_elements import L,J,C,R,Series,Parallel

def get_all_circuit_elements(circuit):
    
    circuit_elements = {}
    capacitors = []
    junctions = []
    inductors = []
    resistors = []
    
    def recursive_function(circuit):
        if type(circuit) == Parallel:
            recursive_function(circuit.component_left)
            recursive_function(circuit.component_right)
        elif type(circuit) == Series:
            recursive_function(circuit.component_left)
            recursive_function(circuit.component_right)
        else:
            circuit_elements[circuit.label] = circuit
            if type(circuit) == L:
                inductors.append(circuit)
            if type(circuit) == J:
                junctions.append(circuit)
            elif type(circuit) == C:
                capacitors.append(circuit)
            elif type(circuit) == R:
                resistors.append(circuit)
    
    recursive_function(circuit)
    return circuit_elements,capacitors,junctions,inductors,resistors


class Bbox(object):
    '''
    Given a circuit, the Bbox allows one to draw the circuit, and compute numerically
    the frequency, dissipation rate and anharmonicity
    of the circuits different modes.
    '''
    def __init__(self, circuit, Q_min = 1.):

        self.circuit = circuit
        self.Q_min = Q_min

        self.all_circuit_elements,\
        self.capacitors,\
        self.junctions,\
        self.inductors,\
        self.resistors = get_all_circuit_elements(circuit)

        self.get_analytical_eigenfrequencies()
        self.circuit_rotated.set_impedence_matrix()
        self.set_analytical_fluxes_anh()
    
    def set_analytical_fluxes_anh(self):
        dY = sp.diff(self.Y,sp.Symbol('w')) # take the derivative
        ImdY = sp.im(dY)
        # Create a function which takes the circuit elements as input and returns eigenfrequencies

        for c in self.all_circuit_elements:   
            comp = self.all_circuit_elements[c]
            flux = comp.flux_wr_ref*sp.sqrt(hbar/sp.Symbol('w')/ImdY)
            voltage = flux*sp.I*sp.Symbol('w')

            comp.flux =\
            sp.utilities.lambdify(
                [elt for elt in self.all_circuit_elements.keys()+["w"]],
                flux,
                "numpy")


            comp.charge =\
            sp.utilities.lambdify(
                [elt for elt in self.all_circuit_elements.keys()+["w"]],
                hbar/2/flux,
                "numpy")


            comp.voltage =\
            sp.utilities.lambdify(
                [elt for elt in self.all_circuit_elements.keys()+["w"]],
                voltage,
                "numpy")


            comp.current =\
            sp.utilities.lambdify(
                [elt for elt in self.all_circuit_elements.keys()+["w"]],
                voltage*comp.admittance(),
                "numpy")

            if type(comp) == J:
                comp.anharmonicity =\
                sp.utilities.lambdify(
                    [elt for elt in self.all_circuit_elements.keys()+["w"]],
                    comp.flux_wr_ref**4*2.*e**2/sp.Symbol(c)/sp.Symbol('w')**2/ImdY**2,
                    "numpy")


    def get_analytical_eigenfrequencies(self):
        # Rotate the circuit such that a junction (or inductor) is in parallel with the rest
        if len(self.junctions)>0:
            self.circuit_rotated = rotate_circuit(self.circuit,self.junctions[0].label) 
            self.reference_LJ = self.junctions[0]
        else:
            self.circuit_rotated = rotate_circuit(self.circuit,self.inductors[0].label) 
            self.reference_LJ = self.inductors[0]

        Y = self.circuit_rotated.admittance()        # Circuit admittance
        Y = sp.together(Y)      # Write as a fraction
        Y_numer = sp.numer(Y)       # Extract the numerator 
        Y_poly = sp.collect(sp.expand(Y_numer),sp.Symbol('w')) # Write numerator as polynomial in omega
        Y_poly_order = sp.polys.polytools.degree(Y_poly,gen = sp.Symbol('w')) # Order of the polynomial
        Y_poly_coeffs = [Y_poly.coeff(sp.Symbol('w'),n) for n in range(Y_poly_order+1)[::-1]] # Get polynomial coefficients
        
        self.Y = Y
        self.N_modes = int(Y_poly_order/2) # Number of modes in the circuit

        # Create a function which takes the circuit elements as input and returns
        # polynomial coeffecients. The root of the polynomial are the 
        # eigenfrequencies and dissipation rates
        self.Y_poly_num = sp.utilities.lambdify(
            [elt for elt in self.all_circuit_elements.keys()],
            Y_poly_coeffs,
            "numpy") # Function which 
        
    def draw(self,
        y_total = 0.6,
        x_total=1.,
        y_C = 0.3,
        x_C = 0.2,
        x_J = 0.2,
        x_LR = 0.7,
        y_LR = 0.3,
        x_par = 0.25,
        fontsize = 15,
        text_position = [0.,-0.35],
        plot = True,
        full_output = False):
        '''
        Draws the circuit
        '''
        line_type = []
                
        def draw(circuit, pos = [0.,0.], which = 'tc'):
            if type(circuit) == Parallel:
                return draw_par(circuit.component_left,circuit.component_right,pos,which)
            elif type(circuit) == Series:
                return draw_ser(circuit.component_left,circuit.component_right,pos,which)
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

            if which == 'l':
                # add side lines

                x+=[np.array([-w/2.,-w/2.-x_par])]
                y+=[np.array([-(h1+h2)/2.,-(h1+h2)/2.])]
                line_type.append('wire')

                x+=[np.array([+w/2.,+w/2.+x_par])]
                y+=[np.array([-(h1+h2)/2.,-(h1+h2)/2.])]
                line_type.append('wire')

                w += x_par*2.

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
        
        element_x,element_y,w,h,xs,ys,element_names = draw(self.circuit,which='t')
        figsize_scaling = 1.5
        fig = plt.figure(figsize = (w*figsize_scaling,h*figsize_scaling))
        ax = fig.add_subplot(111)

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

        element_positions = {}
        for i,el in enumerate(element_names):
            plt.text(element_x[i]+text_position[0],element_y[i]+text_position[1],'$'+el+'$',fontsize=fontsize,ha='center')
            element_positions[el] = [element_x[i],element_y[i]]
        ax.set_axis_off()
        if plot:
            plt.tight_layout()
            plt.show()
        if full_output:
            return element_positions,fig,ax
    
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
        As = np.zeros(len(ws))
        for j in self.junctions:
            As += np.absolute(j.anharmonicity(*args))

        return [ws/2./pi,ks/2./pi,As/h]

    def draw_normal_mode(self,
        mode,unit,
            circuit_parameters,
                fontsize = 15,
                y_arrow=0.25,
                y_arrow_text=0.37):

        def string_to_function(comp,function,circuit_parameters):
            if function == 'flux':
                return comp.flux(**circuit_parameters)/phi_0
            if function == 'charge':
                return comp.charge(**circuit_parameters)/e
            if function == 'voltage':
                return comp.voltage(**circuit_parameters)
            if function == 'current':
                return comp.current(**circuit_parameters)

            
        w,k,A = self.fkA(circuit_parameters)
        element_positions,fig,ax = self.draw(
                plot = False,
                full_output = True,
                y_total = 1,
                x_total=1.,
                y_C = 0.3,
                x_C = 0.2,
                x_J = 0.2,
                x_LR = 0.7,
                y_LR = 0.3,
                fontsize = fontsize,
                text_position = [0.,-0.35])

        circuit_parameters['w'] = (w[mode]+1j*k[mode])*2.*pi
    

        for el,xy in element_positions.items():
            x = xy[0]
            y = xy[1]
            value = string_to_function(self.all_circuit_elements[el],unit,circuit_parameters)
            value_current = string_to_function(self.all_circuit_elements[el],'current',circuit_parameters)
            
            arrow_kwargs = dict(lw = 3,head_width=0.1, head_length=0.1,clip_on=False)
            if np.real(value_current)>0:
                ax.arrow(x-0.25, y+y_arrow,0.5, 0.,fc='red', ec='red',**arrow_kwargs)
            else:
                ax.arrow(x+0.25, y+y_arrow,-0.5, 0.,fc='blue', ec='blue',**arrow_kwargs)
            ax.text(x, y+y_arrow_text,r"%.1e"%np.absolute(value)
                    ,fontsize=fontsize,ha='center')
        plt.show()


if __name__ == '__main__':
    
    circuit = (C('C')|J('L'))|(C('Cc')+(C('Cr')|J('Lr')|R('Rr')))
    b = Bbox(circuit)
    # b.draw()
    kwargs = {}
    kwargs['C'] = 100e-15
    kwargs['L'] = 10e-9
    kwargs['Cc'] = 1e-15
    kwargs['Cr'] = 100e-15
    kwargs['Lr'] = 10e-9
    kwargs['Rr'] = 1e8
    b.draw_normal_mode(0,'current',kwargs)