import os
import core_net as core
from core_net import string_to_component
png_directory = os.path.join(os.path.dirname(__file__),".graphics")

core.pp = {
    "element_width": 1.,
    "element_height": 1.,
    "margin": 0.1,
    "element_height_normal_modes": 1.5,
    "figsize_scaling": 0.5,
    "color": [0.15, 0.15, 0.15],
    "x_fig_margin": 0.2,
    "y_fig_margin": 0.5,
    "C": {
        "gap": 0.2,
        "height": 0.27,
        "lw": 6
    },
    "J": {
        "width": 0.2,
        "lw": 6
    },
    "L": {
        "width": 0.7,
        "height": 0.3,
        "N_points": 150,
        "N_turns": 5,
        "lw": 2
    },
    "R": {
        "width": 0.6,
        "height": 0.35,
        "N_points": 150,
        "N_ridges": 4,
        "lw": 2
    },
    "P": {
        "side_wire_width": 0.25
    },
    "W": {
        "lw": 1
    },
    "label": {
        "fontsize": 10,
        "text_position": 0.35
    },
    "normal_mode_label": {
        "fontsize": 10,
        "y_arrow": 0.26,
        "y_text": 0.37
    },
    "normal_mode_arrow": {
        "logscale": "False",
        "min_width": 0.1,
        "max_width": 0.5,
        "min_lw": 1,
        "max_lw": 3,
        "min_head": 0.07,
        "max_head": 0.071,
        "color_positive": [0.483, 0.622, 0.974],
        "color_negative": [0.931, 0.519, 0.406]
    }
}
# Generate pngs of the different components in their HOVER state
increment = 1
core.pp['W']['lw']+=increment
core.pp['C']['lw']+=increment
core.pp['L']['lw']+=increment
core.pp['R']['lw']+=increment
core.pp['J']['lw']+=increment
try:
    for el in ['R','C','L','J']:
        string_to_component(el,None,None,'').show(save_to = os.path.join(png_directory,'%s_hover.png'%el),plot = False)
except FileNotFoundError:
    os.mkdir(png_directory)
    for el in ['R','C','L','J']:
        string_to_component(el,None,None,'').show(save_to = os.path.join(png_directory,'%s_hover.png'%el),plot = False)

# hover selected state
core.pp['color']=[0.483, 0.622, 0.974]
try:
    for el in ['R','C','L','J']:
        string_to_component(el,None,None,'').show(save_to = os.path.join(png_directory,'%s_hover_selected.png'%el),plot = False)
except FileNotFoundError:
    os.mkdir(png_directory)
    for el in ['R','C','L','J']:
        string_to_component(el,None,None,'').show(save_to = os.path.join(png_directory,'%s_hover_selected.png'%el),plot = False)
# selected state
increment = -1
core.pp['W']['lw']+=increment
core.pp['C']['lw']+=increment
core.pp['L']['lw']+=increment
core.pp['R']['lw']+=increment
core.pp['J']['lw']+=increment
try:
    for el in ['R','C','L','J']:
        string_to_component(el,None,None,'').show(save_to = os.path.join(png_directory,'%s_selected.png'%el),plot = False)
except FileNotFoundError:
    os.mkdir(png_directory)
    for el in ['R','C','L','J']:
        string_to_component(el,None,None,'').show(save_to = os.path.join(png_directory,'%s_selected.png'%el),plot = False)

# rest state
core.pp['color']=[0.15, 0.15, 0.15]
try:
    for el in ['R','C','L','J']:
        string_to_component(el,None,None,'').show(save_to = os.path.join(png_directory,'%s.png'%el),plot = False)
except FileNotFoundError:
    os.mkdir(png_directory)
    for el in ['R','C','L','J']:
        string_to_component(el,None,None,'').show(save_to = os.path.join(png_directory,'%s.png'%el),plot = False)
