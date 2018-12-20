import sympy as sp
import numpy as np
from scipy.constants import e,pi,h,hbar
phi_0 = hbar/2./e
from sympy.core.mul import Mul,Pow,Add
from copy import deepcopy
from json import load
import matplotlib.pyplot as plt
from numbers import Number
from math import floor
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
with open("plotting_parameters.json", "r") as f:
    pp = load(f)

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
        add_vertical_space = False):

        if add_vertical_space:
            element_height_regular = pp['element_height'] 
            pp['element_height'] = pp['element_height_normal_modes']


        element_x,element_y,w,h,xs,ys,element_names,line_type = self.draw(pos = [0.,0.], which = 't')
        fig = plt.figure(figsize = (w*pp["figsize_scaling"],h*pp["figsize_scaling"]))
        ax = fig.add_subplot(111)

        for i in range(len(xs)):
            ax.plot(xs[i],ys[i],color = 'black',lw=pp[line_type[i]]['lw'])

        element_positions = {}
        for i,el in enumerate(element_names):
            plt.text(
                element_x[i]+pp['label']['text_position'][0],
                element_y[i]+pp['label']['text_position'][1],
                el.to_string(),
                fontsize=pp['label']['fontsize']
                ,ha='center')
            element_positions[el] = [element_x[i],element_y[i]]

        if add_vertical_space:
            pp['element_height'] = element_height_regular

        ax.set_axis_off()
        if plot:
            plt.tight_layout()
            plt.show()
        if full_output:
            return element_positions,fig,ax


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
        elif type(circuit.right)==R:
            return self.left
        else:
            return Connection(self.left,self.right)


    #########################
    # Head connection functions
    #########################

    def set_head(self,head = None):
        if self.head is None:
            # The head has not been defined
            if head is None:
                # This connection becomes the head
                self.id_counter = None
                self.Q_min = 1.
                self.Y = None
                self.Y_poly_coeffs_analytical = None
                self.Y_poly_coeffs = None
                self.w_cpx = None
                self.flux_transforms_set = False
                self.dY = None
                self.component_dict = None
                self.inductors = None
                self.capacitors = None
                self.junctions = None
                self.resistors = None

                self.head = self
                self.left.set_head(self)
                self.right.set_head(self)

            else:
                
                self.head = head
                self.left.set_head(head)
                self.right.set_head(head)

        else:
            # The head has already been defined
            pass

    def set_ids(self):
        if self == self.head:
            # First iteration of this recursive function
            if self.head.id_counter is None:
                # We havnt set the ids yet
                self.id_counter = 0

                self._id = self.head.id_counter
                self.head.id_counter += 1
                self.left.set_ids()
                self.right.set_ids()
        else:

            self._id = self.head.id_counter
            self.head.id_counter += 1
            self.left.set_ids()
            self.right.set_ids()

    def set_component_lists(self):

        if self == self.head:
            # First iteration of this recursive function
            if self.component_dict is None:
                # We havnt determined the components yet
                self.component_dict = {}
                self.inductors = []
                self.capacitors = []
                self.junctions = []
                self.resistors = []
                self.no_value_components = []

                self.left.set_component_lists()
                self.right.set_component_lists()  
        else:
            self.left.set_component_lists()
            self.right.set_component_lists()   
    
    def set_circuit_rotated(self):

        self.set_parenthood()
        self.set_head()
        self.set_ids()
        self.set_component_lists()

        if len(self.junctions)>0:
            self.circuit_rotated = rotate(self.junctions[0]) 
        else:
            self.circuit_rotated = rotate(self.inductors[0]) 
    
    def set_Y(self):
        if self.Y is None:
            self.set_circuit_rotated()
            self.Y = self.circuit_rotated.admittance() 

    def set_Y_poly_coeffs(self,**kwargs):
        if self.Y_poly_coeffs_analytical is None:
            self.set_Y()
            Y_together = sp.together(self.Y)
            Y_numer = sp.numer(Y_together)       # Extract the numerator of Y(w)
            Y_poly = sp.collect(sp.expand(Y_numer),sp.Symbol('w')) # Write numerator as polynomial in omega
            Y_poly_order = sp.polys.polytools.degree(Y_poly,gen = sp.Symbol('w')) # Order of the polynomial
            self.Y_poly_coeffs_analytical = [Y_poly.coeff(sp.Symbol('w'),n) for n in range(Y_poly_order+1)[::-1]] # Get polynomial coefficients

        if self.Y_poly_coeffs is None and kwargs is None:
            self.Y_poly_coeffs = [complex(x) for x in self.Y_poly_coeffs_analytical]
        elif kwargs is not None:
            self.Y_poly_coeffs = [complex(x.evalf(subs=kwargs)) for x in self.Y_poly_coeffs_analytical]

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
        self.set_circuit_rotated()
        self.check_kwargs(**kwargs)
        if self.w_cpx is None:
            self.set_Y_poly_coeffs(**kwargs)
            ws_cpx = np.roots(self.Y_poly_coeffs)

            # take only roots with a positive real part (i.e. freq) 
            # and significant Q factors
            relevant_sols = np.argwhere((np.real(ws_cpx)>=0.)&(np.real(ws_cpx)>self.Q_min*np.imag(ws_cpx)))
            ws_cpx=ws_cpx[relevant_sols][:,0]
        
            # Sort solutions with increasing frequency
            order = np.argsort(np.real(ws_cpx))
            self.w_cpx = ws_cpx[order]
        
    def set_dY(self):
        if self.dY is None:
            self.set_Y()
            self.dY = sp.utilities.lambdify(['w']+self.no_value_components,sp.diff(self.Y,sp.Symbol('w')),'numpy')
    
    def compute_all_flux_transformations(self):
        if self.flux_transforms_set == False:
            self.circuit_rotated.set_flux_transforms()
            self.flux_transforms_set = True
    
    def eigenfrequencies(self,**kwargs):
        self.set_w_cpx(**kwargs)
        return np.real(self.w_cpx)/2./pi

    def loss_rates(self,**kwargs):
        self.set_w_cpx(**kwargs)
        return np.imag(self.w_cpx)/2./pi

    def Qs(self,**kwargs):
        self.set_w_cpx(**kwargs)
        return np.real(self.w_cpx)/np.imag(self.w_cpx)
    
    def anharmonicities(self,**kwargs):
        self.set_w_cpx(**kwargs)
        self.set_dY() # the junction flux will be calling dY

        if len(self.junctions) == 0:
            print "There are no junctions and hence no anharmonicity in the circuit"

        elif len(self.junctions) == 1:
            def flux_wr_ref(w,**kwargs):
                return 1.
            self.junctions[0].flux_wr_ref = flux_wr_ref
            return np.absolute(self.junctions[0].anharmonicity(self.w_cpx,**kwargs))/h

        else:
            self.compute_all_flux_transformations()
            return sum([np.absolute(j.anharmonicity(self.w_cpx,**kwargs)) for j in self.junctions])/h

    def show_normal_mode(self,mode,unit='current',**kwargs):

        self.set_w_cpx(**kwargs)
        mode_w = np.real(self.head.w_cpx[mode])
        self.set_dY() # the fluxes will be calling dY
        self.compute_all_flux_transformations()

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


        element_positions,fig,ax = self.show(
        plot = False,
        full_output = True,
        add_vertical_space = True)

        for el,xy in element_positions.items():
            x = xy[0]
            y = xy[1]
            value = string_to_function(el,unit,**kwargs)
            value_current = string_to_function(el,'current',**kwargs)
            
            arrow_kwargs = dict(lw = 3,head_width=0.1, head_length=0.1,clip_on=False)
            if np.real(value_current)>0:
                ax.arrow(x-0.25, y+pp["normal_mode_label"]["y_arrow"],0.5, 0.,fc='red', ec='red',**arrow_kwargs)
            else:
                ax.arrow(x+0.25, y+pp["normal_mode_label"]["y_arrow"],-0.5, 0.,fc='blue', ec='blue',**arrow_kwargs)
            ax.text(x, y+pp["normal_mode_label"]["y_text"],pretty(np.absolute(value),unit)
                    ,fontsize=pp["normal_mode_label"]["fontsize"],ha='center')
        plt.show()

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

    def draw(self,pos,which):
            pos_x1,pos_y1,w1,h1,x1,y1,elt_names1,line_type1 = self.left_o.draw([0.,0.],'l')
            pos_x2,pos_y2,w2,h2,x2,y2,elt_names2,line_type2 = self.right_o.draw([w1,0.],'l')
            w = w1+w2
            h = max(h1,h2)
            x = x1+x2
            y = y1+y2
            line_type = line_type1+line_type2

            if pos == [0.,0.] and which == 't':
                # This is the head connection
                # Add a wire connecting the two sides of the circuit
                y_bottom = -h/2.+pp["label"]["text_position"][1]-pp["margin"]

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

    def draw(self,pos,which):

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
                self.value = a

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

    def draw(self,pos,which):

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
                y_list[i]+=pos[1]-pp['element_height']/2.

        if which == 'l':
            pos_x = pos[0]+pp['element_width']/2.
            pos_y = pos[1]
        elif which == 't':
            pos_x = pos[0]
            pos_y = pos[1]-pp['element_height']/2.

        return [pos_x],[pos_y],pp['element_width'],pp['element_height'],x_list,y_list,[self],line_type


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

    def draw(self,pos,which):
        
        line_type = []
        x = [
            np.array([0.,(pp['element_width']-pp['J']['width'])/2.]),
            np.array([(pp['element_width']+pp['J']['width'])/2.,pp['element_width']]),
            np.array([(pp['element_width']-pp['J']['width'])/2.,(pp['element_width']+pp['J']['width'])/2.]),
            np.array([(pp['element_width']-pp['J']['width'])/2.,(pp['element_width']+pp['J']['width'])/2.]),
            np.array([(pp['element_width']-pp['J']['width'])/2.,(pp['element_width']+pp['J']['width'])/2.])
        ]
        y = [
            np.array([0.,0.]),
            np.array([0.,0.]),
            np.array([0.,0.]),
            np.array([-1.,1.])*pp['J']['width']/2.,
            np.array([1.,-1.])*pp['J']['width']/2.
        ]
        line_type.append('W')
        line_type.append('W')
        line_type.append('W')
        line_type.append('J')
        line_type.append('J')
        
        for i,_ in enumerate(x):
            if which == 'l':
                x[i]+=pos[0]
                y[i]+=pos[1]
                pos_x = pos[0]+pp['element_width']/2.
                pos_y = pos[1]
            elif which == 't':
                x[i]+=pos[0]-pp['element_width']/2.
                y[i]+=pos[1]-pp['element_height']/2.
                pos_x = pos[0]
                pos_y = pos[1]-pp['element_height']/2.

        return [pos_x],[pos_y],pp['element_width'],pp['element_height'],x,y,[self],line_type

class R(Component):
    def __init__(self, arg1 = None, arg2 = None):
        super(R, self).__init__(arg1,arg2)
        self.unit = r'$\Omega$'

    def admittance(self):
        return 1/self.get_value()

    def set_component_lists(self):
        super(R, self).set_component_lists()
        self.head.resistors.append(self)

    def draw(self,pos,which):

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
                y_list[i]+=pos[1]-pp['element_height']/2.

        if which == 'l':
            pos_x = pos[0]+pp['element_width']/2.
            pos_y = pos[1]
        elif which == 't':
            pos_x = pos[0]
            pos_y = pos[1]-pp['element_height']/2.

        return [pos_x],[pos_y],pp['element_width'],pp['element_height'],x_list,y_list,[self],line_type        



class C(Component):
    def __init__(self, arg1 = None, arg2 = None):
        super(C, self).__init__(arg1,arg2)
        self.unit = 'F'
    def admittance(self):
        return sp.I*Mul(sp.Symbol('w'),self.get_value())

    def set_component_lists(self):
        super(C, self).set_component_lists()
        self.head.capacitors.append(self)

    def draw(self,pos,which):
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
            elif which == 't':
                x[i]+=pos[0]-pp['element_width']/2.
                y[i]+=pos[1]-pp['element_height']/2.
                pos_x = pos[0]
                pos_y = pos[1]-pp['element_height']/2.

        return [pos_x],[pos_y],pp['element_width'],pp['element_height'],x,y,[self],line_type
    



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
    rotated_circuit.head = circuit_element.head
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

if __name__ == '__main__':
    circuit = (J(10e-9,'L_{J,1}')|(C(100e-15,'C_J')))+C('C_c')+(C(100e-15)|J(10e-9)|R(1e7))
    # circuit.show()
    circuit.show_normal_mode(0,'charge',C_c = 1e-15)
    circuit.eigenfrequencies(C_c = 1e-15)
    circuit.anharmonicities(C_c = 1e-15)
    # circuit.loss_rates()