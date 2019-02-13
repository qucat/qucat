import sympy as sp
import numpy as np
from scipy.constants import e, pi, h, hbar
from sympy.core.mul import Mul, Pow, Add
from copy import deepcopy
import matplotlib.pyplot as plt
from numbers import Number
from math import floor, factorial
from Qcircuits import gui
import os
from Qcircuits.constants import *
from Qcircuits.utility import pretty_value,\
    check_there_are_no_iterables_in_kwarg, shift, to_string


def string_to_component(s, *arg, **kwarg):
    if s == 'W':
        return W(*arg, **kwarg)
    elif s == 'R':
        return R(*arg, **kwarg)
    elif s == 'L':
        return L(*arg, **kwarg)
    elif s == 'J':
        return J(*arg, **kwarg)
    elif s == 'C':
        return C(*arg, **kwarg)
    elif s == 'G':
        return G(*arg, **kwarg)

plot_parameters = {
    "figsize_scaling": 1,
    "color": light_black,
    "x_fig_margin": 1,
    "y_fig_margin": 0.5,
    "C": {
        "gap": 0.2,
        "height": 0.25,
        "lw": 6
    },
    "J": {
        "width": 0.22,
        "lw": 6
    },
    "L": {
        "width": 0.7,
        "height": 0.25,
        "N_points": 150,
        "N_turns": 5,
        "lw": 2
    },
    "R": {
        "width": 0.6,
        "height": 0.25,
        "N_points": 150,
        "N_ridges": 4,
        "lw": 2
    },
    "P": {
        "side_wire_width": 0.25
    },
    "node": {
        "diameter": 12
    },
    "W": {
        "lw": 1
    },
    "label": {
        "fontsize": 10,
        "text_position_horizontal": [0.,-0.22],
        "text_position_vertical": [0.22,0.]
    }
}

pp = deepcopy(plot_parameters)
scale = 1.25
pp["figsize_scaling"] = scale
pp["C"]["gap"] /= scale
pp["C"]["height"] /= scale
pp["J"]["width"] /= scale
pp["L"]["width"] /= scale
pp["L"]["height"] /= scale
pp["R"]["width"] /= scale
pp["R"]["height"] /= scale
pp["label"]= {
        "fontsize": 10,
        "text_position_horizontal": [0.,-pp["C"]["height"]/2-0.07],
        "text_position_vertical": [pp["C"]["height"]/2+0.05,0.05]
    }
pp["normal_mode_label"]= {
        "fontsize": 10,
        "y_arrow": pp["C"]["height"]/2+0.08,
        "text_position_horizontal": [0.,pp["C"]["height"]/2+0.25],
        "text_position_vertical": [-pp["C"]["height"]/2-0.1,-0.07]
    }
pp["normal_mode_arrow"]= {
        "min_width": 0.1,
        "max_width": 0.4,
        "min_lw": 1,
        "max_lw": 3,
        "min_head": 0.07,
        "max_head": 0.071,
        "color_positive": blue,
        "color_negative": blue
    }
plot_parameters_normal_modes = pp

class _Qcircuit(object):
    """docstring for BBQcircuit"""

    def __init__(self, netlist):
        self.netlist = netlist
        self.network = Network(netlist)
        self.inductors = []
        self.capacitors = []
        self.junctions = []
        self.resistors = []
        self.wires = []
        self.grounds = []
        self.no_value_components = []
        for elt in netlist:
            elt.head = self
            elt.set_component_lists()

        if len(self.junctions) > 0:
            self.ref_elt = self.junctions[0]
        elif len(self.inductors) > 0:
            self.ref_elt = self.inductors[0]
        else:
            raise ValueError(
                "There should be at least one junction or inductor in the circuit")

        self.Q_min = 1.
        self.Y = self.network.admittance(
            self.ref_elt.node_minus, self.ref_elt.node_plus)
        Y_together = sp.together(self.Y)
        Y_numer = sp.numer(Y_together)       # Extract the numerator of Y(w)
        Y_denom = sp.denom(Y_together)       # Extract the numerator of Y(w)
        # Write numerator as polynomial in omega
        Y_numer_poly = sp.collect(sp.expand(Y_numer), sp.Symbol('w'))
        # Write numerator as polynomial in omega
        Y_denom_poly = sp.collect(sp.expand(Y_denom), sp.Symbol('w'))
        self.Y_numer_poly_order = sp.polys.polytools.degree(
            Y_numer_poly, gen=sp.Symbol('w'))  # Order of the polynomial
        self.Y_denom_poly_order = sp.polys.polytools.degree(
            Y_denom_poly, gen=sp.Symbol('w'))  # Order of the polynomial

        self.Y_numer_poly_coeffs_analytical = [Y_numer_poly.coeff(sp.Symbol('w'), n) for n in range(
            self.Y_numer_poly_order+1)[::-1]]  # Get polynomial coefficients
        self.Y_numer_poly_coeffs_analytical = [sp.utilities.lambdify(
            self.no_value_components, c, 'numpy') for c in self.Y_numer_poly_coeffs_analytical]

        self.Y_denom_poly_coeffs_analytical = [Y_denom_poly.coeff(sp.Symbol('w'), n) for n in range(
            self.Y_denom_poly_order+1)[::-1]]  # Get polynomial coefficients
        self.Y_denom_poly_coeffs_analytical = [sp.utilities.lambdify(
            self.no_value_components, c, 'numpy') for c in self.Y_denom_poly_coeffs_analytical]

        self.flux_transformation_dict = {}
        for node in self.network.nodes:
            self.flux_transformation_dict[node] = {}

    def dY(self, w, **kwargs):

        # test if w is an iterable
        try:
            iter(w)
        except TypeError:
            # iterable = False

            # derivative of u/v is (du*v-dv*u)/v^2
            u = sum([complex(a*w**(self.Y_numer_poly_order-n))
                     for n, a in enumerate(self.Y_numer_poly_coeffs(**kwargs))])
            v = sum([complex(a*w**(self.Y_denom_poly_order-n))
                     for n, a in enumerate(self.Y_denom_poly_coeffs(**kwargs))])
            du = sum([complex((self.Y_numer_poly_order-n)*a*w**(self.Y_numer_poly_order-n-1))
                      for n, a in enumerate(self.Y_numer_poly_coeffs(**kwargs))])
            dv = sum([complex((self.Y_denom_poly_order-n)*a*w**(self.Y_denom_poly_order-n-1))
                      for n, a in enumerate(self.Y_denom_poly_coeffs(**kwargs))])

        else:
            # iterable = True

            # derivative of u/v is (du*v-dv*u)/v^2
            u = sum([np.array([complex(a*_w**(self.Y_numer_poly_order-n)) for _w in w])
                     for n, a in enumerate(self.Y_numer_poly_coeffs(**kwargs))])
            v = sum([np.array([complex(a*_w**(self.Y_denom_poly_order-n)) for _w in w])
                     for n, a in enumerate(self.Y_denom_poly_coeffs(**kwargs))])
            du = sum([np.array([complex((self.Y_numer_poly_order-n)*a*_w**(self.Y_numer_poly_order-n-1))
                                for _w in w]) for n, a in enumerate(self.Y_numer_poly_coeffs(**kwargs))])
            dv = sum([np.array([complex((self.Y_denom_poly_order-n)*a*_w**(self.Y_denom_poly_order-n-1))
                                for _w in w]) for n, a in enumerate(self.Y_denom_poly_coeffs(**kwargs))])

        return (du*v-dv*u)/v**2

    def Y_numer_poly_coeffs(self, **kwargs):
        return [complex(coeff(**kwargs)) for coeff in self.Y_numer_poly_coeffs_analytical]

    def Y_denom_poly_coeffs(self, **kwargs):
        return [complex(coeff(**kwargs)) for coeff in self.Y_denom_poly_coeffs_analytical]

    def w_k_A_chi(self, pretty_print=False, **kwargs):

        list_element = None
        list_values = None
        for el, value in kwargs.items():
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
                raise ValueError(
                    "You can only iterate over the value of one element.")

        if pretty_print == True and list_element is not None:
            raise ValueError(
                "Cannot pretty print since $%s$ does not have a unique value" % list_element)

        if list_element is None:

            to_return = self.eigenfrequencies(**kwargs),\
                self.loss_rates(**kwargs),\
                self.anharmonicities(**kwargs),\
                self.kerr(**kwargs)

            if pretty_print:
                N_modes = len(to_return[0])
                table_line = ""
                for i in range(4):
                    table_line += " %7s |"
                table_line += "\n"

                to_print = table_line % (
                    "mode", " freq. ", " diss. ", " anha. ")
                for i, w in enumerate(to_return[0]):
                    to_print += table_line % tuple([str(i)]+[pretty_value(
                        to_return[j][i], use_math=False)+'Hz' for j in range(3)])

                to_print += "\nKerr coefficients\n(diagonal = Kerr, off-diagonal = cross-Kerr)\n"

                table_line = ""
                for i in range(N_modes+1):
                    table_line += " %7s |"
                table_line += "\n"

                to_print += table_line % tuple(['mode'] +
                                               [str(i)+'   ' for i in range(N_modes)])

                for i in range(N_modes):
                    line_elements = [str(i)]
                    for j in range(N_modes):
                        if i >= j:
                            line_elements.append(pretty_value(
                                to_return[3][i][j], use_math=False)+'Hz')
                        else:
                            line_elements.append("")
                    to_print += table_line % tuple(line_elements)
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
            w = np.moveaxis(np.array(w), 0, -1)
            k = np.moveaxis(np.array(k), 0, -1)
            A = np.moveaxis(np.array(A), 0, -1)
            kerr = np.moveaxis(np.array(kerr), 0, -1)
            return w, k, A, kerr

    def check_kwargs(self, **kwargs):
        for key in kwargs:
            if key in self.no_value_components:
                pass
            elif key in [c.label for _, c in self.component_dict.iteritems()]:
                raise ValueError(
                    'The value of %s was already specified when constructing the circuit' % key)
            else:
                raise ValueError(
                    '%s is not the label of a circuit element' % key)

        for label in self.no_value_components:
            try:
                kwargs[label]
            except Exception as e:
                raise ValueError(
                    'The value of %s should be specified with the keyword argument %s=... ' % (label, label))

    def set_w_cpx(self, **kwargs):
        self.check_kwargs(**kwargs)
        ws_cpx = np.roots(self.Y_numer_poly_coeffs(**kwargs))
        
        if len(self.resistors) == 0:
            ws_cpx = np.real(ws_cpx)

        # take only roots with a positive real part (i.e. freq)
        # and significant Q factors
        # Keep solutions with negative imaginary parts, since if there is no resistors in the circuit,
        # there is no garanties on the sign of the resistor...
        relevant_sols = np.argwhere((np.real(ws_cpx) >= 0.) & (
            np.real(ws_cpx) > self.Q_min*np.imag(ws_cpx)))
        ws_cpx = ws_cpx[relevant_sols][:, 0]

        # Sort solutions with increasing frequency
        order = np.argsort(np.real(ws_cpx))
        self.w_cpx = ws_cpx[order]

    def eigenfrequencies(self, **kwargs):
        self.set_w_cpx(**kwargs)
        return np.real(self.w_cpx)/2./pi

    def loss_rates(self, **kwargs):
        self.set_w_cpx(**kwargs)
        return np.imag(self.w_cpx)/2./pi

    def anharmonicities_per_junction(self, pretty_print=False, **kwargs):
        self.set_w_cpx(**kwargs)
        return [j.anharmonicity(self.w_cpx, **kwargs)/h for j in self.junctions]

    def kerr(self, **kwargs):
        As = self.anharmonicities_per_junction(**kwargs)
        N_modes = len(self.w_cpx)
        N_junctions = len(self.junctions)

        Ks = np.zeros((N_modes, N_modes))
        for i in range(N_modes):
            line = []
            for j in range(N_modes):
                for k in range(N_junctions):
                    if i == j:
                        Ks[i, j] += np.absolute(As[k][i])
                    else:
                        Ks[i, j] += 2. * \
                            np.sqrt(np.absolute(As[k][i])
                                    * np.absolute(As[k][j]))
        return Ks

    def anharmonicities(self, **kwargs):
        Ks = self.kerr(**kwargs)
        return [Ks[i, i] for i in range(Ks.shape[0])]

    def hamiltonian(self, 
        modes='all', 
        junc_pot_taylor_exp=4, 
        excitations=6, 
        **kwargs):

        from qutip import destroy, qeye, tensor

        fs = self.eigenfrequencies(**kwargs)
        N_modes = len(fs)
        N_junctions = len(self.junctions)

        if modes == 'all':
            modes = range(N_modes)
        if excitations is not list:
            excitations = [int(excitations) for i in modes]

        H = 0
        phi = [0 for junction in self.junctions]
        qeye_list = [qeye(n) for n in excitations]

        for i, f in enumerate(fs):
            a_list = deepcopy(qeye_list)
            a_list[i] = destroy(excitations[i])
            a = tensor(a_list)
            H += f*a.dag()*a
            phi_0 = hbar/2./e
            for j, junction in enumerate(self.junctions):
                phi[j] += junction.flux(w=f*2.*pi, **kwargs)/phi_0*(a+a.dag())

        for j, junction in enumerate(self.junctions):
            n = 2
            EJ = (hbar/2./e)**2/(junction.get_value(**kwargs)*h)
            while 2*n <= junc_pot_taylor_exp:
                H += (-1)**(n+1)*EJ/factorial(2*n)*phi[j]**(2*n)
                n += 1

        return H


class Qcircuit_NET(_Qcircuit):
    """docstring for Qcircuit"""

    def __init__(self, netlist):
        super(Qcircuit_NET, self).__init__(netlist)


class Qcircuit_GUI(_Qcircuit):
    """docstring for Qcircuit"""

    def __init__(self, filename, edit=True, plot=True, print_network=True):
        
        if edit:
            gui.open_canvas(filename)

        # if file does not exist, also open the gui
        try:
            with open(filename, 'r') as f:
                pass
        except FileNotFoundError as e:
            gui.open_canvas(filename)

        netlist = []
        with open(filename, 'r') as f:
            for el in f:
                el = el.replace('\n', '')
                el = el.split(";")
                if el[3] == '':
                    v = None
                else:
                    v = float(el[3])
                if el[4] == '':
                    l = None
                else:
                    l = el[4]
                netlist.append(
                    string_to_component(el[0], el[1], el[2], v, l))

        super(Qcircuit_GUI, self).__init__(netlist)
        for el in self.netlist:
            el.set_plot_coordinates()

        if plot:
            self.show()

        if print_network:
            for el in [el for el in self.netlist if not isinstance(el,W)]:
                print("%s %d %d %s"%(
                    el.__class__.__name__,
                    min(el.node_minus,el.node_plus),
                    max(el.node_minus,el.node_plus),
                    el.to_string(use_math = False)))
            print('\n')

    def show(self,
             plot=True,
             full_output=False,
             pp=plot_parameters,
             save_to=None,
             **savefig_kwargs):
        
        self.pp = pp

        xs = []
        ys = []
        line_type = []
        for elt in self.netlist:
            x, y, lt = elt.draw()
            xs += x
            ys += y
            line_type += lt

        x_min = min([np.amin(x) for x in xs])
        x_max = max([np.amax(x) for x in xs])
        y_min = min([np.amin(x) for x in ys])
        y_max = max([np.amax(x) for x in ys])

        x_margin = self.pp['x_fig_margin']
        # ensures that any text labels are not cutoff
        y_margin = self.pp['y_fig_margin']
        fig = plt.figure(figsize=(
            ((x_max-x_min)+2.*x_margin)*self.pp["figsize_scaling"],
            ((y_max-y_min)+2.*y_margin)*self.pp["figsize_scaling"]))
        ax = fig.add_subplot(111)
        plt.subplots_adjust(left=0., right=1., top=1., bottom=0.)

        for i, _ in enumerate(xs):
            if line_type[i] == "node":
                ax.scatter(xs[i], ys[i], c=self.pp["color"], s=self.pp['node']['diameter'])
            else:
                ax.plot(xs[i], ys[i], color=self.pp["color"], lw=self.pp[line_type[i]]['lw'])

        for elt in self.netlist:
            elt.draw_label(ax)

        ax.set_axis_off()
        ax.set_xlim(x_min-x_margin, x_max+x_margin)
        ax.set_ylim(y_min-y_margin, y_max+y_margin)
        plt.margins(x=0., y=0.)

        if full_output:
            return fig, ax

        if save_to is not None:
            fig.savefig(save_to, transparent=True, **savefig_kwargs)

        if plot:
            plt.show()

        plt.close()

    def show_normal_mode(self, mode, unit='current',
                         plot=True, save_to=None, **kwargs):

        check_there_are_no_iterables_in_kwarg(**kwargs)
        self.set_w_cpx(**kwargs)
        mode_w = np.real(self.w_cpx[mode])

        def string_to_function(comp, function, **kwargs):
            if function == 'flux':
                phi_0 = hbar/2./e
                return comp.flux(mode_w, **kwargs)/phi_0
            if function == 'charge':
                return comp.charge(mode_w, **kwargs)/e
            if function == 'voltage':
                return comp.voltage(mode_w, **kwargs)
            if function == 'current':
                return comp.current(mode_w, **kwargs)

        def pretty(v, function):
            if function == 'flux':
                return pretty_value(v)+r'$\phi_0$'
            elif function == 'charge':
                return pretty_value(v, use_power_10=True)+r'$e$'
            elif function == 'voltage':
                return pretty_value(v)+'V'
            elif function == 'current':
                return pretty_value(v)+'A'

        fig, ax = self.show(
            plot=False,
            full_output=True,
            pp=plot_parameters_normal_modes)

        # Determine arrow size
        all_values = []
        for el in self.netlist:
            if not isinstance(el,W):
                all_values.append(string_to_function(el, unit, **kwargs))
        all_values = np.absolute(all_values)
        max_value = np.amax(all_values)
        min_value = np.amin(all_values)

        def value_to_01_range(value):
            if max_value == min_value:
                return 1.
            else:
                return (np.absolute(value)-min_value)/(max_value-min_value)

        def arrow_width(value):
            value_01 = value_to_01_range(value)
            ppnm = self.pp['normal_mode_arrow']
            return np.absolute(ppnm['min_width']+value_01*(ppnm['max_width']-ppnm['min_width']))

        def arrow_kwargs(value):
            value_01 = value_to_01_range(value)
            ppnm = self.pp['normal_mode_arrow']
            lw = ppnm['min_lw']+value_01*(ppnm['max_lw']-ppnm['min_lw'])
            head = ppnm['min_head']+value_01 * \
                (ppnm['max_head']-ppnm['min_head'])
            return {'lw': lw,
                    'head_width': head,
                    'head_length': head,
                    'clip_on': False}

        for el in self.netlist:
            if not isinstance(el,W):
                value = string_to_function(el, unit, **kwargs)
                value_current = string_to_function(el, 'current', **kwargs)

                x = el.x_plot_center
                y = el.y_plot_center

                if el.angle%180 == 0.:
                    # Defined for positive arrows
                    x_arrow = x-arrow_width(value)/2.
                    y_arrow = y+self.pp["normal_mode_label"]["y_arrow"]
                    dx_arrow = arrow_width(value)
                    dy_arrow = 0.

                    x_text = x+self.pp["normal_mode_label"]["text_position_horizontal"][0]
                    y_text = y+self.pp["normal_mode_label"]["text_position_horizontal"][1]

                    ha = 'center'
                    va = 'top'

                else:
                    # Defined for positive arrows
                    x_arrow = x-self.pp["normal_mode_label"]["y_arrow"]
                    y_arrow = y-arrow_width(value)/2.
                    dx_arrow = 0.
                    dy_arrow = arrow_width(value)

                    x_text = x+self.pp["normal_mode_label"]["text_position_vertical"][0]
                    y_text = y+self.pp["normal_mode_label"]["text_position_vertical"][1]

                    ha = 'right'
                    va = 'center'

                if np.real(value_current) > 0:
                    ax.arrow(x_arrow, y_arrow, dx_arrow, dy_arrow,
                             fc=self.pp['normal_mode_arrow']['color_positive'],
                             ec=self.pp['normal_mode_arrow']['color_positive'],
                             **arrow_kwargs(value))
                else:
                    ax.arrow(x_arrow+dx_arrow, y_arrow+dy_arrow, -dx_arrow, -dy_arrow,
                             fc=self.pp['normal_mode_arrow']['color_negative'],
                             ec=self.pp['normal_mode_arrow']['color_negative'],
                             **arrow_kwargs(value))

                ax.text(x_text, y_text,
                        pretty(np.absolute(value), unit),
                        fontsize=self.pp["normal_mode_label"]["fontsize"],
                        ha=ha, va=va, style='italic', weight='normal')

        
        w,k,A,chi = self.w_k_A_chi()
        ax.annotate(r'Mode %d, f=%sHz, k=%sHz, A=%sHz'%
            (mode,
            pretty_value(w[mode], use_math=True),
            pretty_value(k[mode], use_math=True),
            pretty_value(A[mode], use_math=True)),
            xy=(0.05, 0.97),
            horizontalalignment='left',
            verticalalignment='center',
            xycoords='axes fraction',
            fontsize=12, 
            weight='bold')

        if plot == True:
            plt.show()
        if save_to is not None:
            fig.savefig(save_to, transparent=True)
        plt.close()


class Network(object):
    """
    The Network class parses network arrays generated by the GUI
    or written manually by the user. 
    It allows the computation of an admittance Y of the network between
    two nodes as well as the transfer function (S_21) of the network
    between two given ports (i.e. 4 nodes).
    Y across an inductor or junction is used to determine the normal-mode 
    frequencies, dissipations and, in conjunction with S_21, 
    the zero-point fluctuations of flux (and other variables) between 
    two arbitrary nodes.

    Parameters
    ----------
    netlist : list
        list of Component objects
    """

    def __init__(self, netlist):

        def merge_chains(chains, i, j):
            '''Merges two chains (two arrays)'''
            to_add = chains[j]
            del chains[j]
            chains[i] = chains[i]+to_add
            return chains

        # ``chains``` is a list of node-lists (node-list = a chain)
        # each chain lists nodes which are connected one to another
        # either through wires or through ground
        chains = []

        # We start by grouping all nodes which are connected by wires 
        # into chains of nodes indexed by a integer which is to become 
        # the new node names  for future calculations
        for el in netlist:

            # Go through all the wires (note: grounds are instances of wires)
            if isinstance(el,W):
                
                # If this element is a ground, call it's negative node '__ground'
                np = el.node_plus
                if type(el) is G:
                    nm = '__ground'
                else:
                    nm = el.node_minus

                # added tells us if both nodes have already 
                # been added to a same chain
                added = False

                # go through all chains to see if it
                # contains the nodes of this wire
                for i, ch in enumerate(chains):

                    # If both node of the wire have already been 
                    # added to a same chain, don't do anyting
                    if (nm in ch) and (np in ch):
                        added = True

                    # If minus node was added to chain 'ch'...
                    elif (nm in ch):
                        for j, ch2 in enumerate(chains):
                            # ...and plus node to chain 'ch2'
                            if np in ch2:
                                # merge the two chains
                                chains = merge_chains(chains, i, j)
                                added = True
                        if added == False:
                            # otherwise add plus node to chain 'ch'
                            ch.append(np)
                            added = True
                    
                    # same check for the plus node
                    elif (np in ch):
                        for j, ch2 in enumerate(chains):
                            if nm in ch2:
                                chains = merge_chains(chains, i, j)
                                added = True
                        if added == False:
                            ch.append(nm)
                            added = True
                
                # if none of the nodes were present in chains, 
                # create a new chain linking the two nodes of the wire
                if added == False:
                    chains.append([nm, np])

        def plot_node_to_new_node(plot_node):
            '''
            Transforms the node ``plot_node``` to a new node.

            Parameters
            ----------
            plot_node:  typically a string or an integer, but could be any
                        hashable object
                        For GUI generated networks, 
                        this is a string 'x,y' that determines the position
                        of the node when plotting it.

            Returns
            -------
            i:          integer, a unique number corresponding to all nodes
                        connected via a wire or ground
        
            '''
            i = 0
            # if plot_node is already in a chain, 
            # return the index of that chain in ``chains``
            for ch in chains:
                if plot_node in ch:
                    return i
                i+=1

            # other wise append ``chains`` and 
            # return the the last index of ``chains``
            chains.append([plot_node])
            return i

        # replace plotting nodes with new nodes
        # for all non-wire elements
        # and make a list of all nodes
        self.nodes = []
        for el in netlist:
            el.node_minus_plot = el.node_minus
            el.node_plus_plot = el.node_plus
            if not isinstance(el,W):
                el.node_minus = plot_node_to_new_node(el.node_minus)
                el.node_plus = plot_node_to_new_node(el.node_plus)
                for n in [el.node_minus, el.node_plus]:
                    if n not in self.nodes:
                        self.nodes.append(n)

        
        # build ``net_dict``, a dictionary such that  
        # ``net_dict[node_A][node_B]`` gives the non-wire circuit 
        # component connecting ``node_A`` and ``node_B``.
        # If ``node_A`` and ``node_B`` are not connected, 
        # calling ``net_dict[node_A][node_B]`` will raise a KeyError
        self.net_dict = {}
        for n in self.nodes:
            self.net_dict[n] = {}
        for el in netlist:
            if not isinstance(el,W):
                self.connect(el, el.node_minus, el.node_plus)

    def connect(self, element, node_minus, node_plus):
            '''
            Modifies the ``net_dict`` variable such that ``node_minus``
            and ``node_plus`` are marked as connected in future calculations.
            ``net_dict`` is a dictionary such that  
            ``net_dict[node_A][node_B]`` gives the non-wire circuit 
            component connecting ``node_A`` and ``node_B``.
            If ``node_A`` and ``node_B`` are not connected, 
            calling ``net_dict[node_A][node_B]`` will raise a KeyError

            Parameters
            ----------
            element:    a ``Circuit`` object which is not a Wire or a Ground
            node_minus: integer
                        negative node of the element
            node_plus: integer
                        positive node of the element        
            '''

        # Connect node minus to node plus
        try:
            # If the nodes are already connected, add this element in parallel 
            self.net_dict[node_minus][node_plus] = self.net_dict[node_minus][node_plus] | element
        except KeyError:
            # Case where the nodes are not connected
            self.net_dict[node_minus][node_plus] = element

        # Connect node plus to node minus
        try:
            self.net_dict[node_plus][node_minus] = self.net_dict[node_plus][node_minus] | element
        except KeyError:
            self.net_dict[node_plus][node_minus] = element

    def remove_node(self, node):

        connections = [x for x in self.net_dict[node].items()]

        # Sum of admittances
        sum_Y = sum([elt.admittance() for _, elt in connections])

        # Prepare mesh
        mesh_to_add = []
        for i, (node_A, elt_A) in enumerate(connections):
            for node_B, elt_B in connections[i+1:]:
                Y = elt_A.admittance()*elt_B.admittance()/sum_Y
                mesh_to_add.append([Admittance(Y), node_A, node_B])
        # Remove star
        for other_node in self.net_dict[node]:
            del self.net_dict[other_node][node]
        del self.net_dict[node]

        # Add mesh
        for mesh_branch in mesh_to_add:
            self.connect(mesh_branch[0], mesh_branch[1], mesh_branch[2])

    def admittance(self, node_minus, node_plus):
        ntr = deepcopy(self) # ntr stands for network to reduce
        for node in self.nodes:
            if node not in [node_minus, node_plus]:
                ntr.remove_node(node)

        Y = ntr.net_dict[node_minus][node_plus].admittance()
        return Y

    def branch_impedance(self, node_1, node_2):
        try:
            return 1./self.net_dict[node_1][node_2].admittance()
        except KeyError:
            return 0.

    def transfer(self, node_left_minus, node_left_plus, node_right_minus, node_right_plus):

        if (node_left_minus in [node_right_plus, node_right_minus])\
            and (node_left_plus in [node_right_plus, node_right_minus]):
            return 1.

        # Reduce network
        ntr = deepcopy(self) # ntr stands for network to reduce
        for node in self.nodes:
            if node not in [node_left_minus, node_left_plus, node_right_minus, node_right_plus]:
                ntr.remove_node(node)

        if (node_left_minus in [node_right_plus, node_right_minus]) or\
            (node_left_plus in [node_right_plus, node_right_minus]):
            if node_left_minus == node_right_minus:
                Z = ntr.branch_impedance(
                    node_left_minus, node_right_minus)
            elif node_left_plus == node_right_plus:
                Z = ntr.branch_impedance(
                    node_left_plus, node_right_plus)

            elif node_left_minus == node_right_plus:
                Z = - \
                    ntr.branch_impedance(
                        node_left_minus, node_right_plus)
            elif node_left_plus == node_right_minus:
                Z = - \
                    ntr.branch_impedance(
                        node_left_plus, node_right_minus)

            # see Pozar
            ABCD = sp.Matrix([[1, Z], [0, 1]])

        else:
            # Compute ABCD of lattice network
            # see https://www.globalspec.com/reference/71734/203279/10-11-lattice-networks
            # Network Analysis & Circuit (By M. Arshad )section 10.11: LATTICE NETWORKS
            Za = ntr.branch_impedance(
                node_left_plus, node_right_plus)
            Zb = ntr.branch_impedance(
                node_left_minus, node_right_plus)
            Zc = ntr.branch_impedance(
                node_left_plus, node_right_minus)
            Zd = ntr.branch_impedance(
                node_left_minus, node_right_minus)
            sum_Z = sum([Za, Zb, Zc, Zd])
            Z11 = (Za+Zb)*(Zd+Zc)/sum_Z
            Z21 = (Zb*Zc-Za*Zd)/sum_Z
            Z22 = (Za+Zc)*(Zd+Zb)/sum_Z

            # see Pozar
            ABCD = sp.Matrix([[
                Z11/Z21,
                Z11*Z22/Z21-Z21], [
                1/Z21,
                Z22/Z21
            ]])

        # Connect missing two elements
        Z_L = ntr.branch_impedance(
            node_left_plus, node_left_minus)
        Z_R = ntr.branch_impedance(
            node_right_plus, node_right_minus)
        ABCD_L = sp.Matrix([[1, Z_L], [0, 1]])
        ABCD_R = sp.Matrix([[1, Z_R], [0, 1]])
        ABCD = ABCD_L*ABCD*ABCD_R

        tr = 1/ABCD[0, 0]
        return tr


class Circuit(object):
    """docstring for Circuit"""

    def __init__(self, node_minus, node_plus):
        self.node_minus = node_minus
        self.node_plus = node_plus
        self.head = None

    def __or__(self, other_circuit):
        return Parallel(self, other_circuit)

    def set_plot_coordinates(self):

        self.x_plot_node_minus = float(self.node_minus_plot.split(',')[0])
        self.x_plot_node_plus = float(self.node_plus_plot.split(',')[0])
        self.y_plot_node_minus = -float(self.node_minus_plot.split(',')[1])
        self.y_plot_node_plus = -float(self.node_plus_plot.split(',')[1])
        self.x_plot_center = (self.x_plot_node_minus +
                              self.x_plot_node_plus)/2.
        self.y_plot_center = (self.y_plot_node_minus +
                              self.y_plot_node_plus)/2.
        if self.x_plot_node_minus == self.x_plot_node_plus:
            # increasing y = SOUTH in tkinter
            if self.y_plot_node_minus < self.y_plot_node_plus:
                self.angle = SOUTH
            else:
                self.angle = NORTH
        else:
            if self.x_plot_node_minus < self.x_plot_node_plus:
                self.angle = WEST
            else:
                self.angle = EAST

    def draw_label(self, ax):
        if self.angle%180. == 0.:
            x = self.x_plot_center+self.head.pp['label']['text_position_horizontal'][0]
            y = self.y_plot_center+self.head.pp['label']['text_position_horizontal'][1]
            ha = 'center'
            va = 'top'

        else:
            x = self.x_plot_center+self.head.pp['label']['text_position_vertical'][0]
            y = self.y_plot_center+self.head.pp['label']['text_position_vertical'][1]
            ha = 'left'
            va = 'center'

        ax.text(x, y,
                to_string(self.unit, self.label, self.value),
                fontsize=self.head.pp['label']['fontsize'],
                ha=ha, va=va)


class Parallel(Circuit):
    """docstring for Connection"""

    def __init__(self, left, right):
        super(Parallel, self).__init__(node_minus=None, node_plus=None)

        # sets the two children circuit elements
        self.left = left
        self.right = right

    def admittance(self):
        return Add(
            self.left.admittance(),
            self.right.admittance())
    
class Component(Circuit):
    """docstring for Component"""

    def __init__(self, node_minus, node_plus, arg1=None, arg2=None):
        super(Component, self).__init__(node_minus, node_plus)
        self.label = None
        self.value = None
        self._flux_wr_ref = None

        if arg1 is None and arg2 is None:
            raise ValueError("Specify either a value or a label")
        for a in [arg1, arg2]:
            if a is None:
                pass
            elif type(a) is str:
                self.label = a
            else:
                self.value = float(a)

    def __hash__(self):
        if self.label is None:
            return hash(str(self.value)+self.unit)
        else:
            if self.value is None:
                return hash(self.label+self.unit)
            else:
                return hash(str(self.value)+self.label+self.unit)

    def get_value(self, **kwargs):
        if self.value is not None:
            return self.value
        elif self.value is None and kwargs is not None:
            if self.label in [k for k in kwargs]:
                return kwargs[self.label]

        return sp.Symbol(self.label)

    def set_component_lists(self):
        if self.value is None and self.label not in ['', ' ', 'None', None]:
            if self.label in self.head.no_value_components:
                raise ValueError(
                    "Two components may not have the same name %s" % self.label)
            else:
                self.head.no_value_components.append(self.label)

    @property
    def flux_wr_ref(self):
        if self._flux_wr_ref is None:
            try:
                tr = self.head.flux_transformation_dict[self.node_minus,
                                                        self.node_plus]
            except KeyError:
                tr = self.head.network.transfer(
                    self.head.ref_elt.node_minus, self.head.ref_elt.node_plus, self.node_minus, self.node_plus)
                self.head.flux_transformation_dict[self.node_minus,
                                                   self.node_plus] = tr
                self.head.flux_transformation_dict[self.node_plus,
                                                   self.node_minus] = -tr

            self._flux_wr_ref = sp.utilities.lambdify(
                ['w']+self.head.no_value_components, tr, "numpy")
        return self._flux_wr_ref

    def flux(self, w, **kwargs):
        ImdY = np.imag(self.head.dY(w, **kwargs))
        return complex(self.flux_wr_ref(w, **kwargs)*sp.sqrt(hbar/w/ImdY))

    def voltage(self, w, **kwargs):
        return complex(self.flux(w, **kwargs)*1j*w)

    def current(self, w, **kwargs):
        kwargs['w'] = w
        Y = self.admittance()
        if isinstance(Y, Number):
            pass
        else:
            Y = Y.evalf(subs=kwargs)
        return complex(self.voltage(**kwargs)*Y)

    def charge(self, w, **kwargs):
        return self.current(w, **kwargs)/w

    def to_string(self, use_math=True, use_unicode=False):
        return to_string(self.unit, self.label, self.value,
                  use_math=use_math, use_unicode=use_unicode)


class W(Component):
    """docstring for Wire"""

    def __init__(self, node_minus, node_plus, arg1='', arg2=None):
        super(W, self).__init__(node_minus, node_plus, arg1='', arg2=None)
        self.unit = None
        self.label = None
        self.value = None

    def to_string(*args, **kwargs):
        return ' '

    def set_component_lists(self):
        super(W, self).set_component_lists()
        self.head.wires.append(self)

    def draw(self):

        x = [np.array([self.x_plot_node_minus, self.x_plot_node_plus]),
        np.array([self.x_plot_node_minus]),
        np.array([self.x_plot_node_plus])]

        y = [np.array([self.y_plot_node_minus, self.y_plot_node_plus]),
        np.array([self.y_plot_node_minus]),
        np.array([self.y_plot_node_plus])]

        line_type = ['W']
        line_type += ['node']
        line_type += ['node']

        return x, y, line_type

class G(W):
    '''
    From a network perspective, a ground element is a wire that connects
    node_plus to a node called '__ground'
    '''
    def __init__(self, node_minus, node_plus, arg1=None, arg2=None):
        super(G, self).__init__(node_minus, node_plus, arg1, arg2)

    def set_component_lists(self):
        super(G, self).set_component_lists()
        self.head.grounds.append(self)

    def draw(self):
        # Defined for EAST
        line_type = []
        x = [
            np.array([0.5, 0.3]),
            np.array([0.3, 0.3]),
            np.array([0.23, 0.23]),
            np.array([0.16, 0.16]),
        ]
        y = [
            np.array([0., 0.]),
            np.array([-1., 1.])*5./30.,
            np.array([-1., 1.])*3./30.,
            np.array([-1., 1.])*1./30.,
        ]
        line_type.append('W')
        line_type.append('W')
        line_type.append('W')
        line_type.append('W')

        if self.angle == WEST:
            return shift(x, self.x_plot_center), shift(y, self.y_plot_center), line_type
        elif self.angle == NORTH:
            return shift(y, self.x_plot_center), shift([-xx for xx in x], self.y_plot_center), line_type
        elif self.angle == EAST:
            return shift([-xx for xx in x], self.x_plot_center), shift(y, self.y_plot_center), line_type
        elif self.angle == SOUTH:
            return shift(y, self.x_plot_center), shift(x, self.y_plot_center), line_type


class R(Component):
    def __init__(self, node_minus, node_plus, arg1=None, arg2=None):
        super(R, self).__init__(node_minus, node_plus, arg1, arg2)
        self.unit = r'$\Omega$'

    def admittance(self):
        return 1/self.get_value()

    def set_component_lists(self):
        super(R, self).set_component_lists()
        self.head.resistors.append(self)

    def draw(self):

        x = np.linspace(-0.25, 0.25 +
                        float(self.head.pp['R']['N_ridges']), self.head.pp['R']['N_points'])
        height = 1.
        period = 1.
        a = height*2.*(-1.+2.*np.mod(np.floor(2.*x/period), 2.))
        b = -height*2.*np.mod(np.floor(2.*x/period), 2.)
        y = (2.*x/period - np.floor(2.*x/period))*a+b+height

        line_type = []
        line_type.append('R')

        # reset leftmost point to 0
        x_min = x[0]
        x -= x_min

        # set width inductor width
        x_max = x[-1]
        x *= self.head.pp['R']['width']/x_max

        # set leftmost point to the length of
        # the side connection wires
        x += (1.-self.head.pp['R']['width'])/2.

        # add side wire connections
        x_min = x[0]
        x_max = x[-1]
        x_list = [x]
        x_list += [np.array([0., x_min])]
        x_list += [np.array([x_max, 1.])]
        line_type.append('W')
        line_type.append('W')

        # center in x
        x_list = shift(x_list-1./2.)

        # set height of inductor
        y *= self.head.pp['R']['height']/2.

        # add side wire connections
        y_list = [y]
        y_list += [np.array([0., 0.])]
        y_list += [np.array([0., 0.])]

        if self.angle%180. == 0.:
            return shift(x_list, self.x_plot_center), shift(y_list, self.y_plot_center), line_type
        if self.angle%180. == 90.:
            return shift(y_list, self.x_plot_center), shift(x_list, self.y_plot_center), line_type

class L(Component):
    def __init__(self, node_minus, node_plus, arg1=None, arg2=None):
        super(L, self).__init__(node_minus, node_plus, arg1, arg2)
        self.unit = 'H'

    def admittance(self):
        return -sp.I*Mul(1/sp.Symbol('w'), 1/self.get_value())

    def set_component_lists(self):
        super(L, self).set_component_lists()
        self.head.inductors.append(self)

    def draw(self):

        x = np.linspace(0.5, float(
            self.head.pp['L']['N_turns']) + 1., self.head.pp['L']['N_points'])
        y = -np.sin(2.*np.pi*x)
        x = np.cos(2.*np.pi*x)+2.*x

        line_type = []
        line_type.append('L')

        # reset leftmost point to 0
        x_min = x[0]
        x -= x_min

        # set width inductor width
        x_max = x[-1]
        x *= self.head.pp['L']['width']/x_max

        # set leftmost point to the length of
        # the side connection wires
        x += (1.-self.head.pp['L']['width'])/2.

        # add side wire connections
        x_min = x[0]
        x_max = x[-1]
        x_list = [x]
        x_list += [np.array([0., x_min])]
        x_list += [np.array([x_max, 1.])]
        line_type.append('W')
        line_type.append('W')

        # center in x
        x_list = shift(x_list, -1./2.)

        # set height of inductor
        y *= self.head.pp['L']['height']/2.

        # add side wire connections
        y_list = [y]
        y_list += [np.array([0., 0.])]
        y_list += [np.array([0., 0.])]

        if self.angle%180. == 0.:
            return shift(x_list, self.x_plot_center), shift(y_list, self.y_plot_center), line_type
        if self.angle%180. == 90.:
            return shift(y_list, self.x_plot_center), shift(x_list, self.y_plot_center), line_type


class J(L):
    def __init__(self, node_minus, node_plus, arg1=None, arg2=None, use_E=False, use_I=False):
        super(J, self).__init__(node_minus, node_plus, arg1, arg2)

        self.use_E = use_E
        self.use_I = use_I
        if self.use_E:
            self.unit = 'Hz'
        elif self.use_I:
            self.unit = 'A'
        else:
            self.unit = 'H'

    def get_value(self, **kwargs):
        value = super(J, self).get_value(**kwargs)
        if (self.use_E == False) and (self.use_I == False):
            return value
        elif (self.use_E == True) and (self.use_I == False):
            L = (hbar/2./e)**2/(value*h)  # E is assumed to be provided in Hz
            return L
        elif (use_E == False) and (use_I == True):
            L = (hbar/2./e)/value
            return L
        else:
            raise ValueError("Cannot set both use_E and use_I to True")

    def set_component_lists(self):
        super(J, self).set_component_lists()
        self.head.junctions.append(self)

    def anharmonicity(self, w, **kwargs):
        ImdY = np.imag(self.head.dY(w, **kwargs))
        return self.flux_wr_ref(w, **kwargs)**4*2.*e**2/self.get_value(**kwargs)/w**2/ImdY**2

    def draw(self):

        line_type = []
        x = [
            np.array([0., 1.]),
            np.array([(1.-self.head.pp['J']['width'])/2.,
                      (1.+self.head.pp['J']['width'])/2.]),
            np.array([(1.-self.head.pp['J']['width'])/2.,
                      (1.+self.head.pp['J']['width'])/2.])
        ]
        y = [
            np.array([0., 0.]),
            np.array([-1., 1.])*self.head.pp['J']['width']/2.,
            np.array([1., -1.])*self.head.pp['J']['width']/2.
        ]
        line_type.append('W')
        line_type.append('J')
        line_type.append('J')

        # center in x and y
        x = shift(x, -1./2.)

        if self.angle%180. == 0.:
            return shift(x, self.x_plot_center), shift(y, self.y_plot_center), line_type
        if self.angle%180. == 90.:
            return shift(y, self.x_plot_center), shift(x, self.y_plot_center), line_type


class R(Component):
    def __init__(self, node_minus, node_plus, arg1=None, arg2=None):
        super(R, self).__init__(node_minus, node_plus, arg1, arg2)
        self.unit = r'$\Omega$'

    def admittance(self):
        return 1/self.get_value()

    def set_component_lists(self):
        super(R, self).set_component_lists()
        self.head.resistors.append(self)

    def draw(self):

        x = np.linspace(-0.25, 0.25 +
                        float(self.head.pp['R']['N_ridges']), self.head.pp['R']['N_points'])
        height = 1.
        period = 1.
        a = height*2.*(-1.+2.*np.mod(np.floor(2.*x/period), 2.))
        b = -height*2.*np.mod(np.floor(2.*x/period), 2.)
        y = (2.*x/period - np.floor(2.*x/period))*a+b+height

        line_type = []
        line_type.append('R')

        # reset leftmost point to 0
        x_min = x[0]
        x -= x_min

        # set width inductor width
        x_max = x[-1]
        x *= self.head.pp['R']['width']/x_max

        # set leftmost point to the length of
        # the side connection wires
        x += (1.-self.head.pp['R']['width'])/2.

        # add side wire connections
        x_min = x[0]
        x_max = x[-1]
        x_list = [x]
        x_list += [np.array([0., x_min])]
        x_list += [np.array([x_max, 1.])]
        line_type.append('W')
        line_type.append('W')

        # center in x
        x_list = shift(x_list, -1./2.)

        # set height of inductor
        y *= self.head.pp['R']['height']/2.

        # add side wire connections
        y_list = [y]
        y_list += [np.array([0., 0.])]
        y_list += [np.array([0., 0.])]

        if self.angle%180. == 0.:
            return shift(x_list, self.x_plot_center), shift(y_list, self.y_plot_center), line_type
        if self.angle%180. == 90.:
            return shift(y_list, self.x_plot_center), shift(x_list, self.y_plot_center), line_type


class C(Component):
    def __init__(self, node_minus, node_plus, arg1=None, arg2=None):
        super(C, self).__init__(node_minus, node_plus, arg1, arg2)
        self.unit = 'F'

    def admittance(self):
        return sp.I*Mul(sp.Symbol('w'), self.get_value())

    def set_component_lists(self):
        super(C, self).set_component_lists()
        self.head.capacitors.append(self)

    def draw(self):
        line_type = []
        x = [
            np.array([0., (1.-self.head.pp['C']['gap'])/2.]),
            np.array([(1.+self.head.pp['C']['gap']) /
                      2., 1.]),
            np.array([(1.-self.head.pp['C']['gap'])/2.,
                      (1.-self.head.pp['C']['gap'])/2.]),
            np.array([(1.+self.head.pp['C']['gap'])/2.,
                      (1.+self.head.pp['C']['gap'])/2.]),
        ]
        y = [
            np.array([0., 0.]),
            np.array([0., 0.]),
            np.array([-self.head.pp['C']['height']/2., self.head.pp['C']['height']/2.]),
            np.array([-self.head.pp['C']['height']/2., self.head.pp['C']['height']/2.]),
        ]
        line_type.append('W')
        line_type.append('W')
        line_type.append('C')
        line_type.append('C')

        # center in x and y
        x = shift(x, -1./2.)

        if self.angle%180. == 0.:
            return shift(x, self.x_plot_center), shift(y, self.y_plot_center), line_type
        if self.angle%180. == 90.:
            return shift(y, self.x_plot_center), shift(x, self.y_plot_center), line_type


class Admittance(Component):
    def __init__(self, Y):
        self.Y = Y

    def admittance(self):
        return self.Y

if __name__ == '__main__':
    c = Qcircuit_GUI('test.txt', edit=True, plot=False, print_network=True)
    c.w_k_A_chi(pretty_print=True)
    c.show_normal_mode(0,unit = 'current')
