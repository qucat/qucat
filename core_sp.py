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

pp={
    "element_width":1.2,
    "element_height":0.6,
    "margin":0.1,
    "element_height_normal_modes":1.0,
    "figsize_scaling":1.,
    "color":[0.15,0.15,0.15],
    "x_fig_margin":0.2,
    "y_fig_margin":0.3,
    "C":{
        "gap":0.2,
        "height":0.27,
        "lw":6
    },
    "J":{
        "width":0.2,
        "lw":6
    },
    "L":{
        "width":0.7,
        "height":0.3,
        "N_points":150,
        "N_turns":5,
        "lw":2
    },
    "R":{
        "width":0.6,
        "height":0.35,
        "N_points":150,
        "N_ridges":4,
        "lw":2
    },
    "P":{
        "side_wire_width":0.25
    },
    "W":{
        "lw":1
    },
    "label":{
        "fontsize":10,
        "text_position":[0.0,-0.35]
    },
    "normal_mode_label":{
        "fontsize":10,
        "y_arrow":0.26,
        "y_text":0.37
    },
    "normal_mode_arrow":{
        "logscale":"False",
        "min_width":0.1,
        "max_width":0.5,
        "min_lw":1,
        "max_lw":3,
        "min_head":0.07,
        "max_head":0.071,
        "color_positive":[0.483, 0.622, 0.974],
        "color_negative":[0.931, 0.519, 0.406]
    }
}

class BBQcircuit(object):
    """docstring for BBQcircuit"""
    def __init__(self, circuit):
        super(BBQcircuit, self).__init__()
        self.circuit = deepcopy(circuit)

        self.circuit.set_head(self)
        self.circuit.set_parenthood()

        self.id_counter = 0
        self.circuit.set_ids()

        self.component_dict = {}
        self.inductors = []
        self.capacitors = []
        self.junctions = []
        self.resistors = []
        self.no_value_components = []
        self.circuit.set_component_lists()

        if len(self.junctions)>0:
            self.circuit_rotated = rotate(self.junctions[0]) 
        else:
            self.circuit_rotated = rotate(self.inductors[0]) 


        self.Q_min = 1.
        self.Y = self.circuit_rotated.admittance() 
        Y_together = sp.together(self.Y)
        Y_numer = sp.numer(Y_together)       # Extract the numerator of Y(w)
        Y_poly = sp.collect(sp.expand(Y_numer),sp.Symbol('w')) # Write numerator as polynomial in omega
        Y_poly_order = sp.polys.polytools.degree(Y_poly,gen = sp.Symbol('w')) # Order of the polynomial
        self.Y_poly_coeffs_analytical = [Y_poly.coeff(sp.Symbol('w'),n) for n in range(Y_poly_order+1)[::-1]] # Get polynomial coefficients
        self.Y_poly_coeffs_analytical = [sp.utilities.lambdify(self.no_value_components,c,'numpy') for c in self.Y_poly_coeffs_analytical]

        self.dY = sp.utilities.lambdify(['w']+self.no_value_components,sp.diff(self.Y,sp.Symbol('w')),'numpy')
        self.circuit_rotated.set_flux_transforms()

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

    def eigenfrequencies(self,**kwargs):
        self.set_w_cpx(**kwargs)
        return np.real(self.w_cpx)
    def loss_rates(self,**kwargs):
        self.set_w_cpx(**kwargs)
        return np.imag(self.w_cpx)
    
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

    def show_normal_mode(self,mode,unit='current',
        plot = True,save_to = None,**kwargs):

        check_there_are_no_iterables_in_kwarg(**kwargs)
        self.set_w_cpx(**kwargs)
        mode_w = np.real(self.w_cpx[mode])

        def string_to_function(comp,function,**kwargs):
            if function == 'flux':
                return comp.flux(mode_w,**kwargs)/phi_0
            if function == 'charge':
                return comp.charge(mode_w,**kwargs)/e
            if function == 'voltage':
                return comp.voltage(mode_w,**kwargs)
            if function == 'current':
                return comp.current(mode_w,**kwargs)

        def pretty(v,function):

            if function == 'flux':
                return pretty_value(v)+r'$\phi_0$'
            if function == 'charge':
                return pretty_value(v,use_power_10 = True)+r'$e$'
            if function == 'voltage':
                return pretty_value(v)+'V'
            if function == 'current':
                return pretty_value(v)+'A'
     
        element_positions,fig,ax = self.circuit.show(
        plot = False,
        full_output = True,
        add_vertical_space = True)

        # Determine arrow size
        all_values = []
        for el in element_positions:
            all_values.append(string_to_function(el,unit,**kwargs))
        all_values = np.absolute(all_values)
        max_value = np.amax(all_values)
        min_value = np.amin(all_values)

        def value_to_01_range(value):
            if pp['normal_mode_arrow']['logscale'] == "True":
                return (np.log10(value)-np.log10(min_value))/(np.log10(max_value)-np.log10(min_value))
            else:
                return (value-min_value)/(max_value-min_value)

        def arrow_width(value):
            value_01 = value_to_01_range(value)
            ppnm = pp['normal_mode_arrow']
            return ppnm['min_width']+value_01*(ppnm['max_width']-ppnm['min_width'])

        def arrow_kwargs(value):
            value_01 = value_to_01_range(value)
            ppnm = pp['normal_mode_arrow']
            lw = ppnm['min_lw']+value_01*(ppnm['max_lw']-ppnm['min_lw'])
            head = ppnm['min_head']+value_01*(ppnm['max_head']-ppnm['min_head'])
            return {'lw':lw,
            'head_width':head,
            'head_length':head,
            'clip_on':False}

        for i,(el,xy) in enumerate(element_positions.items()):
            value = all_values[i]
            value_current = string_to_function(el,'current',**kwargs)
            x = xy[0]
            x_arrow = arrow_width(value)
            y = xy[1]
            y_arrow = y+pp["normal_mode_label"]["y_arrow"]
            
            if np.real(value_current)>0:
                ax.arrow(x-x_arrow/2.,y_arrow ,x_arrow, 0.,
                    fc = pp['normal_mode_arrow']['color_positive'],
                    ec = pp['normal_mode_arrow']['color_positive'],
                    **arrow_kwargs(value))
            else:
                ax.arrow(x+x_arrow/2., y_arrow,-x_arrow, 0.,
                    fc = pp['normal_mode_arrow']['color_negative'],
                    ec = pp['normal_mode_arrow']['color_negative'],
                    **arrow_kwargs(value))
            ax.text(x, y+pp["normal_mode_label"]["y_text"],pretty(value,unit)
                    ,fontsize=pp["normal_mode_label"]["fontsize"],
                    ha='center',style='italic',weight = 'bold')

        if plot == True:
            plt.show()
        if save_to is not None:
            fig.savefig(save_to,transparent = True)
        plt.close()

    def show(*args,**kwargs):
        return self.circuit.show(*args,**kwargs)

class Circuit(object):
    """docstring for Circuit"""
    def __init__(self):
        self._id = None 
        self.parent = None
        self.head = None 

    def __eq__(self, other):
        return self._id == other._id

    def __add__(self, other_circuit):
        return Series(self, other_circuit)

    def __or__(self, other_circuit):
        return Parallel(self, other_circuit)

    def show(self,
        plot = True,
        full_output = False,
        add_vertical_space = False,
        save_to = None):

        if add_vertical_space:
            pp['elt_height'] = pp['element_height_normal_modes']
        else:
            pp['elt_height'] = pp['element_height']



        element_x,element_y,w,h,xs,ys,element_names,line_type = self.draw(pos = [0.,0.], which = 't',is_first_element_to_plot = True)
        x_min = min([np.amin(x) for x in xs])
        x_max = max([np.amax(x) for x in xs])
        y_min = min([np.amin(x) for x in ys])
        y_max = max([np.amax(x) for x in ys])

        x_margin = pp['x_fig_margin']
        y_margin = pp['y_fig_margin'] # ensures that any text labels are not cutoff
        fig = plt.figure(figsize = ((w+2.*x_margin)*pp["figsize_scaling"],(h+2.*y_margin)*pp["figsize_scaling"]))
        ax = fig.add_subplot(111)

        for i in range(len(xs)):
            ax.plot(xs[i],ys[i],color = pp["color"],lw=pp[line_type[i]]['lw'])

        element_positions = {}
        for i,el in enumerate(element_names):
            plt.text(
                element_x[i]+pp['label']['text_position'][0],
                element_y[i]+pp['label']['text_position'][1],
                el.to_string(),
                fontsize=pp['label']['fontsize']
                ,ha='center')
            element_positions[el] = [element_x[i],element_y[i]]

        ax.set_axis_off()
        ax.set_xlim(x_min-x_margin,x_max+x_margin)
        ax.set_ylim(y_min-y_margin,y_max+y_margin)
        plt.margins(x=0.,y=0.)

        if full_output:
            return element_positions,fig,ax

        if save_to is not None:
            fig.savefig(save_to,transparent = True)
        if plot:
            plt.show()
        plt.close()

class Connection(Circuit):
    """docstring for Connection"""
    def __init__(self, left, right):
        super(Connection, self).__init__()

        # sets the two children circuit elements
        self.left = left
        self.right = right

        # o = original
        # stores the children as initially defined
        # (for plotting purposes)
        self.left_o = left
        self.right_o = right

    def set_parenthood(self):
        self.left.parent = self
        self.right.parent = self
        self.left.set_parenthood()
        self.right.set_parenthood()

    def remove_resistances(self):
        if type(self.left)==R:
            return self.right
        elif type(self.right)==R:
            return self.left
        else:
            return Connection(self.left,self.right)

    def set_head(self,head = None):   
        self.head = head
        self.left.set_head(head)
        self.right.set_head(head)


    def set_ids(self):
        self._id = self.head.id_counter
        self.head.id_counter += 1
        self.left.set_ids()
        self.right.set_ids()

    def set_component_lists(self):
        self.left.set_component_lists()
        self.right.set_component_lists()   

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

    def draw(self,pos,which,is_first_element_to_plot = False):
            pos_x1,pos_y1,w1,h1,x1,y1,elt_names1,line_type1 = self.left_o.draw([0.,0.],'l')
            pos_x2,pos_y2,w2,h2,x2,y2,elt_names2,line_type2 = self.right_o.draw([w1,0.],'l')
            w = w1+w2
            h = max(h1,h2)
            x = x1+x2
            y = y1+y2
            line_type = line_type1+line_type2

            if is_first_element_to_plot:
                # Add a wire connecting the two sides of the circuit
                y_extra = (-pp["label"]["text_position"][1]-pp['elt_height']/2.)+pp["margin"]
                y_bottom = -h/2.-y_extra
                h+=y_extra 

                x+=[np.array([0.,0.])]
                y+=[np.array([0.,y_bottom])]
                line_type.append('W')
                x+=[np.array([w,w])]
                y+=[np.array([0.,y_bottom])]
                line_type.append('W')
                x+=[np.array([0.,w])]
                y+=[np.array([y_bottom,y_bottom])]
                line_type.append('W')

            if which == 't':
                return shift(pos_x1+pos_x2,pos[0]-w/2.),shift(pos_y1+pos_y2,pos[1]-h/2.),w,h,shift(x,pos[0]-w/2.),shift(y,pos[1]-h/2.),elt_names1+elt_names2,line_type
            else:
                return shift(pos_x1+pos_x2,pos[0]),shift(pos_y1+pos_y2,pos[1]),w,h,shift(x,pos[0]),shift(y,pos[1]),elt_names1+elt_names2,line_type

class Parallel(Connection):
    """docstring for Parallel"""
    def __init__(self, left,right):
        super(Parallel, self).__init__(left,right)
   
    def admittance(self):
        return Add(self.left.admittance(),
            self.right.admittance())

    def set_flux_transforms(self,ABCD = id2):
        self.right.set_flux_transforms(
            ABCD*sp.Matrix([[1,0],[self.left.admittance(),1]]))
        self.left.set_flux_transforms(
            ABCD*sp.Matrix([[1,0],[self.right.admittance(),1]]))

    def draw(self,pos,which,is_first_element_to_plot = False):

            pos_x1,pos_y1,w1,h1,x1,y1,elt_names1,line_type1 = self.left_o.draw([0.,0.],'t')
            pos_x2,pos_y2,w2,h2,x2,y2,elt_names2,line_type2 = self.right_o.draw([0.,-h1],'t')

            w = max(w1,w2)
            h = h1+h2
            x = x1+x2
            y = y1+y2
            line_type = line_type1+line_type2

            # Connect parallel elements

            x+=[np.array([-w/2.,-w/2.])]
            y+=[np.array([-h1/2.,-h1-h2/2.])]
            line_type.append('W')

            x+=[np.array([+w/2.,+w/2.])]
            y+=[np.array([-h1/2.,-h1-h2/2.])]
            line_type.append('W')

            if np.argmin([w1,w2]) == 0:
                x+=[np.array([-w/2.,-w1/2.])]
                y+=[np.array([-h1/2.,-h1/2.])]
                line_type.append('W')
                x+=[np.array([+w/2.,+w1/2.])]
                y+=[np.array([-h1/2.,-h1/2.])]
                line_type.append('W')
            if np.argmin([w1,w2]) == 1:
                x+=[np.array([-w/2.,-w2/2.])]
                y+=[np.array([-h1-h2/2.,-h1-h2/2.])]
                line_type.append('W')
                x+=[np.array([+w/2.,+w2/2.])]
                y+=[np.array([-h1-h2/2.,-h1-h2/2.])]
                line_type.append('W')

            if which == 'l':

                # add side lines

                x+=[np.array([-w/2.,-w/2.-pp['P']["side_wire_width"]])]
                y+=[np.array([-(h1+h2)/2.,-(h1+h2)/2.])]
                line_type.append('W')

                x+=[np.array([+w/2.,+w/2.+pp['P']["side_wire_width"]])]
                y+=[np.array([-(h1+h2)/2.,-(h1+h2)/2.])]
                line_type.append('W')

                w += pp['P']["side_wire_width"]*2.

            # Move into position

            if which == 'l':
                return shift(pos_x1+pos_x2,pos[0]+w/2.),shift(pos_y1+pos_y2,pos[1]+h/2.),w,h,shift(x,pos[0]+w/2.),shift(y,pos[1]+h/2.),elt_names1+elt_names2,line_type
            else:
                return shift(pos_x1+pos_x2,pos[0]),shift(pos_y1+pos_y2,pos[1]),w,h,shift(x,pos[0]),shift(y,pos[1]),elt_names1+elt_names2,line_type


class Component(Circuit):
    """docstring for Component"""
    def __init__(self, arg1 = None, arg2 = None):
        super(Component, self).__init__()
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

    def remove_resistances(self):
        return self

    def set_parenthood(self):
        pass

    def set_head(self,head = None):
        self.head = head

    def set_ids(self):
        self._id = self.head.id_counter
        self.head.id_counter += 1

    def set_component_lists(self):
        self.head.component_dict[self._id] = self
        if self.value is None:
            if self.label in self.head.no_value_components:
                raise ValueError("Two components may not have the same name %s"%self.label)
            else:
                self.head.no_value_components.append(self.label)

    def set_flux_transforms(self,ABCD = id2):
        ABCD = ABCD*sp.Matrix([[1,0],[self.admittance(),1]])
        self.flux_wr_ref = sp.utilities.lambdify(['w']+self.head.no_value_components,1/ABCD[0,0],"numpy")

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
    def __init__(self, arg1 = None, arg2 = None):
        super(L,self).__init__(arg1,arg2)
        self.unit = 'H'
    def admittance(self):
        return -sp.I*Mul(1/sp.Symbol('w'),1/self.get_value())

    def set_component_lists(self):
        super(L, self).set_component_lists()
        self.head.inductors.append(self)

    def draw(self,pos,which,is_first_element_to_plot = False):

        x = np.linspace(0.5,float(pp['L']['N_turns']) +1. ,pp['L']['N_points'])
        y = -np.sin(2.*np.pi*x)
        x = np.cos(2.*np.pi*x)+2.*x

        line_type = []
        line_type.append('L')
            
        # reset leftmost point to 0
        x_min = x[0]
        x -= x_min

        # set width inductor width
        x_max = x[-1]
        x *= pp['L']['width']/x_max

        # set leftmost point to the length of 
        # the side connection wires
        x +=(pp['element_width']-pp['L']['width'])/2.   

        # add side wire connections
        x_min = x[0]
        x_max = x[-1]
        x_list = [x]
        x_list += [np.array([0.,x_min])]
        x_list += [np.array([x_max,pp['element_width']])]
        line_type.append('W')
        line_type.append('W')

        # Shift into position in x 
        for i,x in enumerate(x_list):
            if which == 'l':
                x_list[i]+=pos[0]
            elif which == 't':
                x_list[i]+=pos[0]-pp['element_width']/2.


        # set height of inductor
        y *=pp['L']['height']/2.

        # add side wire connections
        y_list = [y]
        y_list += [np.array([0.,0.])]
        y_list += [np.array([0.,0.])]

        for i,y in enumerate(y_list):
            if which == 'l':
                y_list[i]+=pos[1]
            elif which == 't':
                y_list[i]+=pos[1]-pp['elt_height']/2.

        if which == 'l':
            pos_x = pos[0]+pp['element_width']/2.
            pos_y = pos[1]
            height = pp['L']['height']
        elif which == 't':
            pos_x = pos[0]
            pos_y = pos[1]-pp['elt_height']/2.
            height = pp['elt_height']

        return [pos_x],[pos_y],pp['element_width'],height,x_list,y_list,[self],line_type

class J(L):
    def __init__(self, arg1 = None, arg2 = None,use_E=False,use_I=False):
        super(J,self).__init__(arg1,arg2)

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

    def draw(self,pos,which,is_first_element_to_plot = False):
        
        line_type = []
        x = [
            np.array([0.,pp['element_width']]),
            np.array([(pp['element_width']-pp['J']['width'])/2.,(pp['element_width']+pp['J']['width'])/2.]),
            np.array([(pp['element_width']-pp['J']['width'])/2.,(pp['element_width']+pp['J']['width'])/2.])
        ]
        y = [
            np.array([0.,0.]),
            np.array([-1.,1.])*pp['J']['width']/2.,
            np.array([1.,-1.])*pp['J']['width']/2.
        ]
        line_type.append('W')
        line_type.append('J')
        line_type.append('J')
        
        for i,_ in enumerate(x):
            if which == 'l':
                x[i]+=pos[0]
                y[i]+=pos[1]
                pos_x = pos[0]+pp['element_width']/2.
                pos_y = pos[1]
                height = pp['J']['width']
            elif which == 't':
                x[i]+=pos[0]-pp['element_width']/2.
                y[i]+=pos[1]-pp['elt_height']/2.
                pos_x = pos[0]
                pos_y = pos[1]-pp['elt_height']/2.
                height = pp['elt_height']

        return [pos_x],[pos_y],pp['element_width'],height,x,y,[self],line_type

class R(Component):
    def __init__(self, arg1 = None, arg2 = None):
        super(R, self).__init__(arg1,arg2)
        self.unit = r'$\Omega$'

    def admittance(self):
        return 1/self.get_value()

    def set_component_lists(self):
        super(R, self).set_component_lists()
        self.head.resistors.append(self)

    def draw(self,pos,which,is_first_element_to_plot = False):

        x = np.linspace(-0.25,0.25+float(pp['R']['N_ridges']),pp['R']['N_points'])
        height = 1.
        period = 1.
        a = height*2.*(-1.+2.*np.mod(np.floor(2.*x/period),2.))
        b = -height*2.*np.mod(np.floor(2.*x/period),2.)
        y = (2.*x/period - np.floor(2.*x/period))*a+b+height

        line_type = []
        line_type.append('R')
            
        # reset leftmost point to 0
        x_min = x[0]
        x -= x_min

        # set width inductor width
        x_max = x[-1]
        x *= pp['R']['width']/x_max

        # set leftmost point to the length of 
        # the side connection wires
        x +=(pp['element_width']-pp['R']['width'])/2.   

        # add side wire connections
        x_min = x[0]
        x_max = x[-1]
        x_list = [x]
        x_list += [np.array([0.,x_min])]
        x_list += [np.array([x_max,pp['element_width']])]
        line_type.append('W')
        line_type.append('W')

        # Shift into position in x 
        for i,x in enumerate(x_list):
            if which == 'l':
                x_list[i]+=pos[0]
            elif which == 't':
                x_list[i]+=pos[0]-pp['element_width']/2.


        # set height of inductor
        y *=pp['R']['height']/2.

        # add side wire connections
        y_list = [y]
        y_list += [np.array([0.,0.])]
        y_list += [np.array([0.,0.])]

        for i,y in enumerate(y_list):
            if which == 'l':
                y_list[i]+=pos[1]
            elif which == 't':
                y_list[i]+=pos[1]-pp['elt_height']/2.

        if which == 'l':
            pos_x = pos[0]+pp['element_width']/2.
            pos_y = pos[1]
            height = pp['R']['height']
        elif which == 't':
            pos_x = pos[0]
            pos_y = pos[1]-pp['elt_height']/2.
            height = pp['elt_height']

        return [pos_x],[pos_y],pp['element_width'],height,x_list,y_list,[self],line_type        

class C(Component):
    def __init__(self, arg1 = None, arg2 = None):
        super(C, self).__init__(arg1,arg2)
        self.unit = 'F'
    def admittance(self):
        return sp.I*Mul(sp.Symbol('w'),self.get_value())

    def set_component_lists(self):
        super(C, self).set_component_lists()
        self.head.capacitors.append(self)

    def draw(self,pos,which,is_first_element_to_plot = False):
        line_type = []
        x = [
            np.array([0.,(pp['element_width']-pp['C']['gap'])/2.]),
            np.array([(pp['element_width']+pp['C']['gap'])/2.,pp['element_width']]),
            np.array([(pp['element_width']-pp['C']['gap'])/2.,(pp['element_width']-pp['C']['gap'])/2.]),
            np.array([(pp['element_width']+pp['C']['gap'])/2.,(pp['element_width']+pp['C']['gap'])/2.]),
        ]
        y = [
            np.array([0.,0.]),
            np.array([0.,0.]),
            np.array([-pp['C']['height']/2.,pp['C']['height']/2.]),
            np.array([-pp['C']['height']/2.,pp['C']['height']/2.]),
        ]
        line_type.append('W')
        line_type.append('W')
        line_type.append('C')
        line_type.append('C')
        
        for i,_ in enumerate(x):
            if which == 'l':
                x[i]+=pos[0]
                y[i]+=pos[1]
                pos_x = pos[0]+pp['element_width']/2.
                pos_y = pos[1]
                height = pp['C']['height']
            elif which == 't':
                x[i]+=pos[0]-pp['element_width']/2.
                y[i]+=pos[1]-pp['elt_height']/2.
                pos_x = pos[0]
                pos_y = pos[1]-pp['elt_height']/2.
                height = pp['elt_height']

        return [pos_x],[pos_y],pp['element_width'],height,x,y,[self],line_type

    


def rotate(circuit_element):
    
    # Check if the tree is already correctly rotated
    if circuit_element.parent.parent is None:
        return circuit_element.parent
    
    def recursive_function(elt):
        parent_connection = elt.parent
        if parent_connection.parent is None:
            if parent_connection.left==elt:
                return parent_connection.right
            elif parent_connection.right == elt:
                return parent_connection.left
        else:
            if parent_connection.left == elt:
                parent_connection.left = recursive_function(parent_connection)
            elif parent_connection.right == elt:
                parent_connection.right = recursive_function(parent_connection)
            return parent_connection
        
    rotated_circuit = circuit_element|recursive_function(circuit_element)
    rotated_circuit._id = 0
    rotated_circuit.parent = None
    rotated_circuit.set_parenthood()

    return rotated_circuit

def shift(to_shift,shift):
    for i,_ in enumerate(to_shift):
        to_shift[i]+= shift
    return to_shift

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
    qubit = C(100e-15)|J(1e-9)
    resonator = C(100e-15)|L(10e-9)|R(1e6)
    cQED_circuit = BBQcircuit(qubit + C(1e-15) + resonator)
    cQED_circuit.show_normal_mode(0)
    cQED_circuit.w_k_A_chi()