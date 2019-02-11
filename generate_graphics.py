import os
import Qcircuits.core_net as core
from Qcircuits.core_net import string_to_component
from Qcircuits.constants import *
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
png_directory = os.path.join(os.path.dirname(__file__),".graphics")
try:
    os.mkdir(png_directory)
except FileExistsError:
    pass
dpi = 300

def generate_icon(comp,hover = False, selected = False):
    pp = {
        "element_width": 1.,
        "element_height": 1.,
        "margin": 0.,
        "figsize_scaling": 1,
        "color": [0.15, 0.15, 0.15],
        "x_fig_margin": 0.,
        "y_fig_margin": 0.25,
        "C": {
            "gap": 0.2,
            "height": 0.25,
            "lw": 6
        },
        "J": {
            "width": 0.25,
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
        "W": {
            "lw": 2
        }
    }
    core.pp = pp

    comp.node_minus_plot = '0,0'
    comp.node_plus_plot = '1,0'
    comp.set_plot_coordinates()

    xs, ys, line_type = comp.draw()


    fig = plt.figure(figsize=(1,0.5))
    ax = fig.add_subplot(111)
    ax.set_axis_off()
    plt.margins(x=0., y=0.)
    ax.set_ylim(-0.25, 0.25)
    ax.set_xlim(0., 1.)
    plt.subplots_adjust(left=0., right=1., top=1., bottom=0.)

    rect_args = [(0.1,-0.22),0.8,0.44]
    rect_kwargs = {'linewidth':2,'facecolor':'none'}

    if hover and selected:
        increment = 1.5
        core.pp['W']['lw']+=increment
        core.pp['C']['lw']+=increment
        core.pp['L']['lw']+=increment
        core.pp['R']['lw']+=increment
        core.pp['J']['lw']+=increment
        state_string = '_hover_selected'
        ax.add_patch(Rectangle(*rect_args,edgecolor=blue,**rect_kwargs))
    elif hover:
        increment = 1.5
        core.pp['W']['lw']+=increment
        core.pp['C']['lw']+=increment
        core.pp['L']['lw']+=increment
        core.pp['R']['lw']+=increment
        core.pp['J']['lw']+=increment
        state_string = '_hover'
        # ax.add_patch(Rectangle(*rect_args,edgecolor=lighter_blue,**rect_kwargs))
    elif selected:
        state_string = '_selected'
        ax.add_patch(Rectangle(*rect_args,edgecolor=light_blue,**rect_kwargs))
    else:
        state_string = ''

    for i in range(len(xs)):
        ax.plot(xs[i], ys[i], color=pp["color"], lw=pp[line_type[i]]['lw'])
    
    fig.savefig(os.path.join(png_directory, comp.__class__.__name__+state_string+'.png'), transparent=True, dpi = dpi)
    plt.close()

for el in ['R','C','L','J','G']:
    generate_icon(string_to_component(el,None,None,''))
    generate_icon(string_to_component(el,None,None,''),hover = True)
    generate_icon(string_to_component(el,None,None,''),selected = True)
    generate_icon(string_to_component(el,None,None,''),hover = True,selected = True)
