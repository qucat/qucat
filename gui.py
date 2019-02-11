try:
    # Tkinter for Python 2.xx
    import Tkinter as tk
    from tkFont import Font
    from Tkinter import tkMessageBox as messagebox
except ImportError:
    # Tkinter for Python 3.xx
    import tkinter as tk
    from tkinter.font import Font
    from tkinter import messagebox
from PIL import Image, ImageTk
from tkinter import ttk
import numpy as np
import os
from Qcircuits.utility import to_string
from Qcircuits.constants import *
from copy import deepcopy

png_directory = os.path.join(os.path.dirname(__file__), ".graphics")
node_dot_radius = 1./30.
lw = 1./50.
lw_hover = 2.*lw
lw_select_hover = 5.*lw
lw_select = 3.*lw


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


class AutoScrollbar(ttk.Scrollbar):
    """ A scrollbar that hides itself if it's not needed. """

    def set(self, lo, hi):
        if float(lo) <= 0.0 and float(hi) >= 1.0:
            self.grid_remove()
        else:
            self.grid()
            ttk.Scrollbar.set(self, lo, hi)

    def pack(self, **kw):
        raise tk.TclError(
            'Cannot use pack with the widget ' + self.__class__.__name__)

    def place(self, **kw):
        raise tk.TclError(
            'Cannot use place with the widget ' + self.__class__.__name__)


class SnappingCanvas(tk.Canvas):
    def __init__(self, master, grid_unit, netlist_file, **kw):
        """ Initialize the ImageFrame """

        self.frame = ttk.Frame()
        self.frame.grid()  # place Canvas widget on the grid
        self.frame.grid(sticky='nswe')  # make frame container sticky
        self.frame.rowconfigure(0, weight=1)  # make canvas expandable
        self.frame.columnconfigure(0, weight=1)

        menu_label_template = "{:<6}"
        # TODO make this pretty?
        menu_font = Font(family="Courier New", size=9, weight='normal')

        self.menubar = tk.Menu(self.frame)
        menu = tk.Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(
            label=menu_label_template.format("File"), menu=menu)

        label_template = "{:<15}{:>6}"
        menu.add_command(label=label_template.format(
            "Save", "Ctrl+S"), command=self.save, font=menu_font)
        menu.add_command(label=label_template.format(
            "Exit", ""), command=master.destroy, font=menu_font)

        menu = tk.Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(
            label=menu_label_template.format("Edit"), menu=menu)

        menu.add_command(label=label_template.format(
            "Undo", "Ctrl+Z"), command=self.ctrl_z, font=menu_font)
        menu.add_command(label=label_template.format(
            "Redo", "Ctrl+Y"), command=self.ctrl_y, font=menu_font)
        menu.add_separator()
        menu.add_command(label=label_template.format(
            "Cut", "Ctrl+X"), command=self.cut_selection, font=menu_font)
        menu.add_command(label=label_template.format(
            "Copy", "Ctrl+C"), command=self.copy_selection, font=menu_font)
        menu.add_command(label=label_template.format(
            "Paste", "Ctrl+V"), command=self.paste, font=menu_font)
        menu.add_separator()
        menu.add_command(label=label_template.format(
            "Select all", "Ctrl+A"), command=self.select_all, font=menu_font)
        menu.add_separator()
        menu.add_command(label=label_template.format(
            "Delete", "Del"), command=self.delete_selection, font=menu_font)
        menu.add_command(label=label_template.format(
            "Delete all", ""), command=self.delete_all, font=menu_font)
        master.config(menu=self.menubar)

        menu = tk.Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(
            label=menu_label_template.format("Insert"), menu=menu)

        menu.add_command(label=label_template.format("Wire", "W"), command=(
            lambda: self.event_generate('w')), font=menu_font)
        menu.add_command(label=label_template.format("Junction", "J"), command=(
            lambda: self.event_generate('j')), font=menu_font)
        menu.add_command(label=label_template.format("Inductor", "L"), command=(
            lambda: self.event_generate('l')), font=menu_font)
        menu.add_command(label=label_template.format("Capacitor", "C"), command=(
            lambda: self.event_generate('c')), font=menu_font)
        menu.add_command(label=label_template.format("Resistor", "R"), command=(
            lambda: self.event_generate('r')), font=menu_font)
        master.config(menu=self.menubar)

        # Vertical and horizontal scrollbars for canvas
        hbar = AutoScrollbar(self.frame, orient='horizontal')
        vbar = AutoScrollbar(self.frame, orient='vertical')
        hbar.grid(row=1, column=0, sticky='we')
        vbar.grid(row=0, column=1, sticky='ns')

        tk.Canvas.__init__(self, self.frame, bd=0, highlightthickness=0,
                           xscrollcommand=hbar.set, yscrollcommand=vbar.set, confine=False, bg="white")

        self.grid(row=0, column=0, sticky='nswe')

        hbar.configure(command=self.scroll_x)  # bind scrollbars to the canvas
        vbar.configure(command=self.scroll_y)
        self.update()  # Waits for the canvas to pop up before asking for its window size below

        self.canvas_center = [
            self.canvasx(self.winfo_width()/2.),
            self.canvasy(self.winfo_height()/2.)]
        self.netlist_file = netlist_file
        self.grid_unit = int(grid_unit)

        self.focus_set()
        self.bind('r', lambda event: R(self, event))
        self.bind('l', lambda event: L(self, event))
        self.bind('c', lambda event: C(self, event))
        self.bind('j', lambda event: J(self, event))
        self.bind('w', lambda event: W(self, event))
        self.bind('g', lambda event: G(self, event))
        self.bind('<Control-s>', self.save)
        self.bind("<Configure>", self.on_resize)
        self.bind('<Delete>', self.delete_selection)
        self.bind('<Control-c>', self.copy_selection)
        self.bind('<Control-x>', self.cut_selection)
        self.bind('<Control-v>', self.paste)
        self.bind('<Control-a>', self.select_all)
        self.bind('<Control-y>', self.ctrl_y)
        self.bind('<Control-z>', self.ctrl_z)
        # zoom for Windows and MacOS, but not Linux
        self.bind('<Control-MouseWheel>', self.wheel)
        # zoom for Linux, wheel scroll down
        self.bind('<Control-Button-5>',   self.wheel)
        # zoom for Linux, wheel scroll up
        self.bind('<Control-Button-4>',   self.wheel)
        # zoom for Windows and MacOS, but not Linux
        self.bind('<Shift-MouseWheel>', self.scroll_x_wheel)
        # zoom for Linux, wheel scroll down
        self.bind('<Shift-Button-5>',   self.scroll_x_wheel)
        # zoom for Linux, wheel scroll up
        self.bind('<Shift-Button-4>',   self.scroll_x_wheel)
        # zoom for Windows and MacOS, but not Linux
        self.bind('<MouseWheel>', self.scroll_y_wheel)
        # zoom for Linux, wheel scroll down
        self.bind('<Button-5>',   self.scroll_y_wheel)
        # zoom for Linux, wheel scroll up
        self.bind('<Button-4>',   self.scroll_y_wheel)

        self.elements = []
        self.copied_elements = []
        self.history = []
        self.history_location = -1
        self.track_changes = False
        try:
            with open(netlist_file, 'r') as f:
                netlist_file_string = [line for line in f]
        except FileNotFoundError:
            netlist_file_string = []
            with open(netlist_file, 'w') as f:
                pass
        self.load_netlist(netlist_file_string)
        self.draw_grid()
        self.configure_scrollregion()
        self.track_changes = True
        self.save()

        # Keep track of what the user knows what to do
        self.used_arrows = False

    def cut_selection(self, event=None):
        self.track_changes = False
        self.copied_elements = [deepcopy(el)
                                for el in self.elements if el.selected]
        self.track_changes = True
        self.delete_selection()

    def paste(self, event=None):
        if len(self.copied_elements) > 0:
            self.deselect_all()

            # smallest x and y of copied elements, in canvas units
            x_min = min([el.x_minus for el in self.copied_elements] +
                        [el.x_plus for el in self.copied_elements])
            y_min = min([el.y_minus for el in self.copied_elements] +
                        [el.y_plus for el in self.copied_elements])
            x_min, y_min = self.grid_to_canvas([x_min, y_min])

            # shift to apply, in canvas units
            dx = self.canvasx(event.x)-x_min
            dy = self.canvasy(event.y)-y_min

            for el in self.copied_elements:
                el.create()
                el.adapt_to_grid_unit()
                el.force_select()
                el.move(dx, dy)
                el.add_or_replace_label()

            self.bind("<Motion>", el.on_motion)
            self.bind("<ButtonPress-1>", el.release_motion_paste)

    def copy_selection(self, event=None):
        self.track_changes = False
        self.copied_elements = [deepcopy(el)
                                for el in self.elements if el.selected]
        self.track_changes = True

    def scroll_y_wheel(self, event):
        if event.num == 5 or event.delta < 0:
            direction = 1
        if event.num == 4 or event.delta > 0:
            direction = -1
        self.yview_scroll(direction, tk.UNITS)
        self.configure_scrollregion()
        self.draw_grid(event)

    def scroll_x_wheel(self, event):
        if event.num == 5 or event.delta < 0:
            direction = 1
        if event.num == 4 or event.delta > 0:
            direction = -1
        self.xview_scroll(direction, tk.UNITS)
        self.configure_scrollregion()
        self.draw_grid(event)

    def wheel(self, event):
        old_grid_unit = self.grid_unit

        scaling = 1.08
        try:
            if abs(event.delta) > 120:
                # Fast scrolling case on Windows
                scaling = 1.15
        except:
            pass
        smallest_grid_unit = 35
        largest_grid_unit = 100

        # Respond to Linux (event.num) or Windows (event.delta) wheel event
        if event.num == 5 or event.delta < 0:  # scroll down, smaller
            new_grid_unit = int(self.grid_unit/scaling)
            if new_grid_unit == old_grid_unit:
                new_grid_unit -= 1
        elif event.num == 4 or event.delta > 0:  # scroll up, bigger
            new_grid_unit = int(self.grid_unit*scaling)
            if new_grid_unit == old_grid_unit:
                new_grid_unit += 1

        if smallest_grid_unit > new_grid_unit:
            self.message("Can't zoom out more")
        elif new_grid_unit > largest_grid_unit:
            self.message("Can't zoom in more")
        else:
            grid_pos_old = self.canvas_to_grid(
                [self.canvasx(event.x), self.canvasy(event.y)])
            self.grid_unit = new_grid_unit
            canvas_pos_old = self.grid_to_canvas(grid_pos_old)
            canvas_center_shift = [self.canvasx(
                event.x)-canvas_pos_old[0], self.canvasy(event.y)-canvas_pos_old[1]]
            self.canvas_center = [self.canvas_center[0]+canvas_center_shift[0],
                                  self.canvas_center[1]+canvas_center_shift[1]]
            for el in self.elements:
                el.adapt_to_grid_unit()

            self.draw_grid(event)
            self.configure_scrollregion()

    def scroll_x(self, *args, **kwargs):
        """ Scroll canvas horizontally and redraw the image """
        self.xview(*args)  # scroll vertically
        self.draw_grid()

    def scroll_y(self, *args, **kwargs):
        """ Scroll canvas vertically and redraw the image """
        self.yview(*args)  # scroll vertically
        self.draw_grid()

    def ctrl_z(self, event=None):
        # print('CTRL-Z')
        # for net in self.history:
        #     print(net)
        # print("location was: %d" % self.history_location)
        if self.history_location > 0:
            self.track_changes = False
            self.history_location -= 1
            self.load_netlist(self.history[self.history_location].split('\n'))
            self.save()
            self.track_changes = True
        else:
            self.message('Nothing to undo')
        # print("location is: %d" % self.history_location)

    def ctrl_y(self, event=None):
        # print('CTRL-Y')
        # for net in self.history:
        #     print(net)
        # print("location was: %d" % self.history_location)
        if 0 <= self.history_location < len(self.history)-1:
            self.track_changes = False
            self.history_location += 1
            self.load_netlist(self.history[self.history_location].split('\n'))
            self.track_changes = True
        else:
            self.message('Nothing to redo')
        # print("location is: %d" % self.history_location)

    def load_netlist(self, lines):
        self.delete_all(track_changes=False)
        for el in lines:
            el = el.replace('\n', '')
            el = el.split(";")
            string_to_component(el[0], self, auto_place=el)

    def on_resize(self, event=None):
        self.configure_scrollregion()
        self.draw_grid(event)

    def configure_scrollregion(self):

        box_canvas = [self.canvasx(0),  # get visible area of the canvas
                      self.canvasy(0),
                      self.canvasx(self.winfo_width()),
                      self.canvasy(self.winfo_height())]

        if len(self.elements) > 0:
            xs = [el.x_minus for el in self.elements] + \
                [el.x_plus for el in self.elements]
            ys = [el.y_minus for el in self.elements] + \
                [el.y_plus for el in self.elements]
            box_elements = self.grid_to_canvas(
                [min(xs)-1, min(ys)-1])+self.grid_to_canvas([max(xs)+1, max(ys)+1])

            self.configure(
                scrollregion=[min(box_elements[0], box_canvas[0]), min(box_elements[1], box_canvas[1]),
                              max(box_elements[2], box_canvas[2]), max(box_elements[3], box_canvas[3])])
        else:
            self.configure(
                scrollregion=box_canvas)

    def draw_grid(self, event=None):

        # Draw the grid
        self.delete("grid")
        dx = 1
        dy = 1
        box_canvas = (self.canvasx(0),  # get visible area of the canvas
                      self.canvasy(0),
                      self.canvasx(self.winfo_width()),
                      self.canvasy(self.winfo_height()))

        self.background = self.create_rectangle(
            *box_canvas, fill='white', outline='', tags='grid')

        grid_x = np.arange(
            self.canvas_center[0], box_canvas[2], self.grid_unit).tolist()
        grid_x += np.arange(self.canvas_center[0]-self.grid_unit,
                            box_canvas[0], -self.grid_unit).tolist()
        grid_y = np.arange(
            self.canvas_center[1], box_canvas[3], self.grid_unit).tolist()
        grid_y += np.arange(self.canvas_center[1]-self.grid_unit,
                            box_canvas[1], -self.grid_unit).tolist()

        for x in grid_x:
            for y in grid_y:
                self.create_line(x-dx, y, x+2*dx, y, tags='grid')
                self.create_line(x, y-dy, x, y+2*dy, tags='grid')
        self.tag_lower('grid')
        self.tag_bind('grid', '<ButtonPress-1>', self.start_selection_field)
        self.tag_bind('grid', "<B1-Motion>", self.expand_selection_field)
        self.tag_bind('grid', "<ButtonRelease-1>", self.end_selection_field)
        self.tag_bind('grid', "<Button-3>", self.right_click)

    def right_click(self, event):
        self.deselect_all()
        self.bind("<ButtonRelease-3>", self.open_right_click_menu)

    def open_right_click_menu(self, event):
        menu = tk.Menu(self, tearoff=0)
        menu.add_command(label="Paste", command=self.paste)
        menu.tk_popup(event.x_root, event.y_root, 0)
        self.bind("<ButtonRelease-3>", lambda event: None)

    def start_selection_field(self, event):
        self.deselect_all()
        self.selection_rectangle_x_start = self.canvasx(event.x)
        self.selection_rectangle_y_start = self.canvasy(event.y)
        self.selection_rectangle = self.create_rectangle(
            self.canvasx(event.x), self.canvasy(event.y), self.canvasx(event.x), self.canvasy(event.y), dash=(3, 5))

    def expand_selection_field(self, event):
        self.deselect_all()
        self.coords(self.selection_rectangle,
                    min(self.canvasx(event.x), self.selection_rectangle_x_start),
                    min(self.canvasy(event.y), self.selection_rectangle_y_start),
                    max(self.canvasx(event.x), self.selection_rectangle_x_start),
                    max(self.canvasy(event.y), self.selection_rectangle_y_start))
        for el in self.elements:
            el.box_select(*self.coords(self.selection_rectangle))

    def end_selection_field(self, event):
        self.delete(self.selection_rectangle)

    def deselect_all(self, event=None):
        for el in self.elements:
            el.deselect()

    def select_all(self, event=None):
        for el in self.elements:
            el.force_select()

    def delete_selection(self, event=None, track_changes=None):
        was_tracking_changes = self.track_changes
        self.track_changes = False
        to_delete = [el for el in self.elements if el.selected]
        for el in to_delete:
            el.delete()
        if track_changes is None:  # Just follow "was_tracking_changes"
            self.track_changes = was_tracking_changes
            self.save()
        elif track_changes is False:
            self.track_changes = was_tracking_changes
        elif track_changes is True:
            self.track_changes = True
            self.save()
            self.track_changes = track_changes
        self.configure_scrollregion()

    def delete_all(self, event=None, track_changes=None):
        self.select_all()
        self.delete_selection(event, track_changes)

    def save(self, event=None):

        netlist_string = ""
        for el in self.elements:
            v, l = el.prop
            if v is None:
                v = ''
            else:
                v = "%e" % v

            if l is None:
                l = ''

            netlist_string += ("%s;%s;%s;%s;%s\n" % (
                type(el).__name__,
                el.grid_to_node_string(el.x_minus, el.y_minus),
                el.grid_to_node_string(el.x_plus, el.y_plus),
                v, l))

        with open(self.netlist_file, 'w') as f:
            f.write(netlist_string)

        if self.track_changes:
            del self.history[self.history_location+1:]
            self.history.append(netlist_string)
            self.history_location += 1
            # print('ADDED TO HISTORY')
            # for net in self.history:
            #     print(net)

        self.message("Saving...")

    def message(self, text, t = 0.3):
        saved_message = self.create_text(
            self.canvasx(5), self.canvasy(2), text=text, anchor=tk.NW,
            font=Font(family='Helvetica', size=8, weight='normal'))
        self.after(int(1000*t), lambda: self.delete(saved_message))

    def create_circle(self, x, y, r):
        # Everything defined in canvas units
        x0 = x - r
        y0 = y - r
        x1 = x + r
        y1 = y + r
        return self.create_oval(x0, y0, x1, y1, fill='black')

    def update_circle(self, circle, x, y, r):
        # Everything defined in canvas units
        x0 = x - r
        y0 = y - r
        x1 = x + r
        y1 = y + r
        self.coords(circle, x0, y0, x1, y1)

    def grid_to_canvas(self, pos):
        # pos = [x ,y] (in grid units)
        return [self.canvas_center[0]+self.grid_unit*pos[0],
                self.canvas_center[1]+self.grid_unit*pos[1]]

    def canvas_to_grid(self, pos):
        # pos = [x ,y] (in canvas units)
        return [(pos[0]-self.canvas_center[0])/self.grid_unit,
                (pos[1]-self.canvas_center[1])/self.grid_unit]

class TwoNodeElement(object):
    def __init__(self, canvas, event=None, auto_place=None):
        self.canvas = canvas
        self.was_moved = False
        self.was_rotated = False
        self.hover = False
        self.selected = False
        self.dot_minus = None
        self.dot_plus = None

        if auto_place is None and event is not None:
            self.manual_place(event)

        elif auto_place is not None and event is None:
            v = auto_place[3]
            l = auto_place[4]
            if l == '':
                l = None

            if v == '':
                v = None
            else:
                v = float(v)

            self.prop = [v, l]

            self.x_minus, self.y_minus = self.node_string_to_grid(
                auto_place[1])
            self.x_plus, self.y_plus = self.node_string_to_grid(
                auto_place[2])
            self.auto_place(auto_place)

    @property
    def pos(self):
        return []

    @pos.setter
    def pos(self, pos):
        pass

    @property
    def prop(self):
        return []

    @prop.setter
    def prop(self, prop):
        pass

    def __deepcopy__(self, memo):
        cls = self.__class__
        newone = cls.__new__(cls)
        newone.__init__(self.canvas)
        memo[id(self)] = newone
        newone.pos = deepcopy(self.pos)
        newone.prop = deepcopy(self.prop)
        return newone

    def set_nodes(self):
        pass

    def grid_to_node_string(self, x, y):
        gu = self.canvas.grid_unit
        return "%d,%d" % (int(x), int(y))

    def node_string_to_grid(self, node):
        xy = node.split(',')
        x = int(xy[0])
        y = int(xy[1])
        return x, y

    def grid_to_canvas(self, pos):
        return self.canvas.grid_to_canvas(pos)

    def deselect(self):
        pass

    def force_select(self):
        pass

    def on_click(self, event, shift_control=False):

        if self.selected is False and shift_control is False:
            self.canvas.deselect_all()

        self.canvas.bind(
            '<Left>', lambda event: self.on_updownleftright(event, angle=WEST))
        self.canvas.bind(
            '<Right>', lambda event: self.on_updownleftright(event, angle=EAST))
        self.canvas.bind(
            '<Up>', lambda event: self.on_updownleftright(event, angle=NORTH))
        self.canvas.bind(
            '<Down>', lambda event: self.on_updownleftright(event, angle=SOUTH))

    def on_motion(self, event):
        x, y = self.get_center_pos()
        dx = self.canvas.canvasx(event.x) - x
        dy = self.canvas.canvasy(event.y) - y
        for el in self.canvas.elements:
            if el.selected or el == self:
                el.move(dx, dy)
        self.was_moved = True

    def release_motion(self, event, shift_control=False):
        N_selected = 0
        self.canvas.track_changes = False


        for el in self.canvas.elements:
            if el.selected or el == self:
                N_selected += 1
                el.snap_to_grid()
                el.add_or_replace_label()

        if self.was_rotated:
            self.set_nodes()
            self.was_rotated = False

        self.canvas.track_changes = True
        self.canvas.save()
        self.canvas.bind('<Left>', lambda event: None)
        self.canvas.bind('<Right>', lambda event: None)
        self.canvas.bind('<Up>', lambda event: None)
        self.canvas.bind('<Down>', lambda event: None)

        if self.was_moved:
            self.force_select()
            self.was_moved = False
        elif shift_control:
            self.ctrl_shift_select()
        else:
            self.select()

    def release_motion_paste(self, event):
        self.release_motion(event)
        self.canvas.bind("<Motion>", lambda event: None)
        self.canvas.bind("<ButtonPress-1>", lambda event: None)

    def hover_enter(self, event):
        self.hover = True
        self.update_graphic()

        self.canvas.tag_bind(self.binding_object, "<Button-1>", self.on_click)
        self.canvas.tag_bind(self.binding_object, "<Shift-Button-1>",
                             lambda event: self.on_click(event, shift_control=True))
        self.canvas.tag_bind(self.binding_object, "<Control-Button-1>",
                             lambda event: self.on_click(event, shift_control=True))
        self.canvas.tag_bind(self.binding_object,
                             "<B1-Motion>", self.on_motion)
        self.canvas.tag_bind(
            self.binding_object, "<ButtonRelease-1>", self.release_motion)
        self.canvas.tag_bind(
            self.binding_object, '<Double-Button-1>', self.double_click)
        self.canvas.tag_bind(self.binding_object, "<Shift-ButtonRelease-1>",
                             lambda event: self.release_motion(event, shift_control=True))
        self.canvas.tag_bind(self.binding_object, "<Control-ButtonRelease-1>",
                             lambda event: self.release_motion(event, shift_control=True))

    def set_allstate_bindings(self):
        self.canvas.tag_bind(self.binding_object, "<Enter>", self.hover_enter)
        self.canvas.tag_bind(self.binding_object, "<Leave>", self.hover_leave)
        self.canvas.tag_bind(self.binding_object,
                             "<Button-3>", self.right_click)

    def hover_leave(self, event):
        self.hover = False
        self.update_graphic()

    def right_click(self, event):
        self.canvas.bind("<ButtonRelease-3>", self.open_right_click_menu)

    def select(self):
        self.canvas.deselect_all()
        if self.selected is False:
            self.selected = True
            self.update_graphic()

    def ctrl_shift_select(self):
        if self.selected is False:
            self.selected = True
            self.update_graphic()
        elif self.selected is True:
            self.deselect()

    def force_select(self):
        self.selected = True
        self.update_graphic()

    def deselect(self):
        self.selected = False
        self.update_graphic()
        
    def add_nodes(self, to = 'all wires'):

        # Check if one of the nodes of this component/wire 
        # is on a wire
        if to == 'all wires':
            to_nodify = []
            for el in self.canvas.elements:
                if type(el) == W and el != self:
                    to_nodify.append(el)
        else:
            to_nodify = to
        
        for w in to_nodify:
            if w.x_minus == w.x_plus == self.x_minus:
                # vertical wire, minus node
                if self.y_minus in range(min(w.y_minus,w.y_plus)+1,
                                    max(w.y_minus,w.y_plus)):
                    x = w.x_minus
                    ym = w.y_minus
                    yp = w.y_plus
                    w.delete()
                    W(self.canvas,auto_place=['W','%d,%d'%(x,ym),
                        '%d,%d'%(x,self.y_minus),'',''])
                    W(self.canvas,auto_place=['W','%d,%d'%(x,yp),
                        '%d,%d'%(x,self.y_minus),'',''])
                    return True

            elif w.y_minus == w.y_plus == self.y_minus:
                # horizontal wire, minus node
                if self.x_minus in range(min(w.x_minus,w.x_plus)+1,
                                    max(w.x_minus,w.x_plus)):
                    xm = w.x_minus
                    xp = w.x_plus
                    y = w.y_minus
                    w.delete()
                    W(self.canvas,auto_place=['W','%d,%d'%(xm,y),
                        '%d,%d'%(self.x_minus,y),'',''])
                    W(self.canvas,auto_place=['W','%d,%d'%(xp,y),
                        '%d,%d'%(self.x_minus,y),'',''])
                    return True
            elif w.x_minus == w.x_plus == self.x_plus:
                # vertical wire, positive node
                if self.y_plus in range(min(w.y_minus,w.y_plus)+1,
                                    max(w.y_minus,w.y_plus)):
                    x = w.x_minus
                    ym = w.y_minus
                    yp = w.y_plus
                    w.delete()
                    W(self.canvas,auto_place=['W','%d,%d'%(x,ym),
                        '%d,%d'%(x,self.y_plus),'',''])
                    W(self.canvas,auto_place=['W','%d,%d'%(x,yp),
                        '%d,%d'%(x,self.y_plus),'',''])
                    return True

            elif w.y_minus == w.y_plus == self.y_plus:
                # horizontal wire, positive node
                if self.x_plus in range(min(w.x_minus,w.x_plus)+1,
                                    max(w.x_minus,w.x_plus)):
                    xm = w.x_minus
                    xp = w.x_plus
                    y = w.y_minus
                    w.delete()
                    W(self.canvas,auto_place=['W','%d,%d'%(xm,y),
                        '%d,%d'%(self.x_plus,y),'',''])
                    W(self.canvas,auto_place=['W','%d,%d'%(xp,y),
                        '%d,%d'%(self.x_plus,y),'',''])
                    return True
        return False

    def add_or_replace_node_dots(self):
        gu = self.canvas.grid_unit
        canvas_coords_minus = self.grid_to_canvas([self.x_minus,self.y_minus])
        canvas_coords_plus = self.grid_to_canvas([self.x_plus,self.y_plus])


        if self.dot_minus is None:
            self.dot_minus = self.canvas.create_circle(
                *canvas_coords_minus, gu*node_dot_radius)
        else:
            self.canvas.update_circle(self.dot_minus,
                *canvas_coords_minus, gu*node_dot_radius)


        if self.dot_plus is None:
            self.dot_plus = self.canvas.create_circle(
                *canvas_coords_plus, gu*node_dot_radius)
        else:
            self.canvas.update_circle(self.dot_plus,
                *canvas_coords_plus, gu*node_dot_radius)


class W(TwoNodeElement):

    def __init__(self, canvas, event=None, auto_place=None):
        self.x_minus = None
        self.y_minus = None
        self.x_plus = None
        self.y_plus = None
        super(W, self).__init__(canvas, event, auto_place)

    @property
    def binding_object(self):
        return self.line

    @property
    def pos(self):
        return [self.x_minus, self.y_minus, self.x_plus, self.y_plus]

    @pos.setter
    def pos(self, pos):
        if pos != self.pos:
            self.x_minus = pos[0]
            self.y_minus = pos[1]
            self.x_plus = pos[2]
            self.y_plus = pos[3]
            self.canvas.save()

    @property
    def prop(self):
        return [None, None]

    @prop.setter
    def prop(self, prop):
        pass
    def add_nodes(self,to = 'all wires'):

        all_other_elements = [el for el in self.canvas.elements if el != self]
        for el in all_other_elements:
            added_a_node = TwoNodeElement.add_nodes(el,to = [self])
            if added_a_node:
                return True

        added_a_node = TwoNodeElement.add_nodes(self,to)
        if added_a_node:
            return True
        else:
            return False

    def box_select(self, x0, y0, x1, y1):
        xs = [x0, x1]
        ys = [y0, y1]
        xm,ym = self.canvas.grid_to_canvas([self.x_minus,self.y_minus])
        xp,yp = self.canvas.grid_to_canvas([self.x_plus,self.y_plus])

        if min(xs) <= xm <= max(xs) and min(ys) <= ym <= max(ys)\
                and min(xs) <= xp <= max(xs) and min(ys) <= yp <= max(ys):
            self.force_select()

    def get_center_pos(self):
        # returns center in canvas units
        xm, ym, xp, yp = self.canvas.coords(self.line)
        return [(xm+xp)/2., (ym+yp)/2.]

    def open_right_click_menu(self, event):
        menu = tk.Menu(self.canvas, tearoff=0)
        menu.add_command(label="Delete", command=self.canvas.delete_selection)
        menu.add_command(label="Copy", command=self.canvas.copy_selection)
        menu.add_command(label="Cut", command=self.canvas.cut_selection)
        menu.tk_popup(event.x_root, event.y_root, 0)
        self.canvas.bind("<ButtonRelease-3>", lambda event: None)

    def double_click(self, event=None):
        pass

    def snap_to_grid(self):
        # in canvas coordinates
        xm, ym, xp, yp = self.canvas.coords(self.line)
        # snapped to grid units
        xm, ym, xp, yp = [round(p) for p in self.canvas.canvas_to_grid(
            [xm, ym])+self.canvas.canvas_to_grid([xp, yp])]
        self.pos = [xm, ym, xp, yp]

        # back to canvas units:
        xm, ym, xp, yp = self.canvas.grid_to_canvas(
            [xm, ym])+self.canvas.grid_to_canvas([xp, yp])
        self.canvas.coords(self.line, xm, ym, xp, yp)
        self.canvas.update_circle(
            self.dot_minus, xm, ym, self.canvas.grid_unit*node_dot_radius)
        self.canvas.update_circle(
            self.dot_plus, xp, yp, self.canvas.grid_unit*node_dot_radius)

        # Add node if a node of this line is on another line
        self.add_nodes()

    def init_minus_snap_to_grid(self, event):
        gu = float(self.canvas.grid_unit)
        x0, y0 = self.canvas.canvas_center
        return int(round(float(self.canvas.canvasx(event.x)-x0)/gu)),\
            int(round(float(self.canvas.canvasy(event.y)-y0)/gu))

    def init_plus_snap_to_grid(self, event):

        xm, ym = self.grid_to_canvas([self.x_minus, self.y_minus])
        xp = self.canvas.canvasx(event.x)
        yp = self.canvas.canvasy(event.y)
        gu = float(self.canvas.grid_unit)
        x0, y0 = self.canvas.canvas_center

        if abs(xm-xp) > abs(ym-yp):
            # Horizontal line
            self.x_plus = int(round(float(self.canvas.canvasx(event.x)-x0)/gu))
            self.y_plus = self.y_minus
        else:
            # Vertical line
            self.x_plus = self.x_minus
            self.y_plus = int(round(float(self.canvas.canvasy(event.y)-y0)/gu))
        
    def end_line(self, event):
        self.canvas.delete("temp")
        self.canvas.bind("<Button-1>", lambda event: None)
        self.canvas.bind('<Motion>', lambda event: None)
        self.canvas.config(cursor='arrow')

        self.init_plus_snap_to_grid(event)
        self.create()
        self.track_changes = False
        self.add_nodes()
        self.track_changes = True
        self.canvas.save()

    def move(self, dx, dy):
        '''
        Input given in canvas units
        '''
        self.canvas.move(self.line, dx, dy)
        self.canvas.move(self.dot_minus, dx, dy)
        self.canvas.move(self.dot_plus, dx, dy)

    def add_or_replace_label(self):
        pass

    def update_graphic(self):
        if self.selected and self.hover:
            self.canvas.itemconfig(self.line, fill=blue, width=lw_select_hover*self.canvas.grid_unit)
        elif self.selected:
            self.canvas.itemconfig(self.line, fill=light_blue, width=lw_select*self.canvas.grid_unit)
        elif self.hover:
            self.canvas.itemconfig(self.line, fill=light_black, width=lw_hover*self.canvas.grid_unit)
        else:
            self.canvas.itemconfig(self.line, fill=light_black, width=lw*self.canvas.grid_unit)

    def on_updownleftright(self, event=None, angle = None):
        pass

    def manual_place(self, event):
        self.canvas.config(cursor='plus')
        self.canvas.bind("<Button-1>", self.start_line)
        self.canvas.bind("<Escape>", self.abort_creation)

    def auto_place(self, auto_place_info):
        self.create()
        self.add_nodes()

    def start_line(self, event):
        self.x_minus, self.y_minus = self.init_minus_snap_to_grid(event)

        
        gu = self.canvas.grid_unit
        self.dot_minus = self.canvas.create_circle(
            gu*self.x_minus, gu*self.y_minus, gu*node_dot_radius)

        self.canvas.bind("<Motion>", self.show_line)
        self.canvas.bind("<Button-1>", self.end_line)

    def abort_creation(self, event=None):
        self.canvas.bind("<Button-1>", lambda event: None)
        self.canvas.bind("<Escape>", lambda event: None)
        self.canvas.bind("<Motion>", lambda event: None)
        self.canvas.bind("<Button-1>", lambda event: None)
        self.canvas.delete('temp')
        self.canvas.config(cursor='arrow')
        del self

    def delete(self, event=None):
        self.canvas.elements.remove(self)
        self.canvas.delete(self.line)
        self.canvas.delete(self.dot_minus)
        self.canvas.delete(self.dot_plus)
        self.canvas.save()
        del self

    def create(self):
        canvas_coords_minus = self.grid_to_canvas(self.pos[:2])
        canvas_coords_plus = self.grid_to_canvas(self.pos[2:])
        self.line = self.canvas.create_line(
            *(canvas_coords_minus+canvas_coords_plus),
            width=lw*self.canvas.grid_unit,
            fill = light_black)
        self.add_or_replace_node_dots()
        self.canvas.elements.append(self)
        self.set_allstate_bindings()

    def adapt_to_grid_unit(self):
        canvas_coords_minus = self.grid_to_canvas(self.pos[:2])
        canvas_coords_plus = self.grid_to_canvas(self.pos[2:])
        self.canvas.coords(
            self.line, *(canvas_coords_minus+canvas_coords_plus))
        self.update_graphic()
        self.add_or_replace_node_dots()


    def show_line(self, event):
        self.canvas.delete("temp")

        xm, ym = self.grid_to_canvas([self.x_minus, self.y_minus])
        xp = self.canvas.canvasx(event.x)
        yp = self.canvas.canvasy(event.y)

        if abs(xm-xp) > abs(ym-yp):
            # Horizontal line
            self.canvas.create_line(xm, ym, xp, ym, tags='temp',
            width=lw*self.canvas.grid_unit,
            fill = light_black)
        else:
            # Vertical line
            self.canvas.create_line(xm, ym, xm, yp, tags='temp',
            width=lw*self.canvas.grid_unit,
            fill = light_black)



class Component(TwoNodeElement):
    def __init__(self, canvas, event=None, auto_place=None):
        self.image = None
        self._value = None
        self._label = None
        self.text = None
        self._x_center = None
        self._y_center = None
        self._angle = None
        super(Component, self).__init__(canvas, event, auto_place)

    def get_center_pos(self):
        # returns center in canvas units
        return self.canvas.coords(self.image)

    @property
    def binding_object(self):
        return self.image

    @property
    def pos(self):
        return [self._x_center, self._y_center, self._angle]

    @pos.setter
    def pos(self, pos):
        # Defined in grid units, _x/y_center should be n+(0. or 0.5) with n an integer
        if pos != self.pos:
            self._x_center = pos[0]
            self._y_center = pos[1]
            self._angle = pos[2]
            self.set_nodes()
            self.canvas.save()

    @property
    def prop(self):
        return [self._value, self._label]

    @prop.setter
    def prop(self, prop):
        if prop != self.prop:
            self._value = prop[0]
            self._label = prop[1]
            self.canvas.save()

    def set_nodes(self):

        # positive y points south
        pos = self.pos
        if self._angle == SOUTH:
            self.x_minus = pos[0]
            self.y_minus = pos[1]-0.5
            self.x_plus = pos[0]
            self.y_plus = pos[1]+0.5
        elif self._angle == NORTH:
            self.x_minus = pos[0]
            self.y_minus = pos[1]+0.5
            self.x_plus = pos[0]
            self.y_plus = pos[1]-0.5

        elif self._angle == EAST:
            self.x_minus = pos[0]-0.5
            self.y_minus = pos[1]
            self.x_plus = pos[0]+0.5
            self.y_plus = pos[1]
        elif self._angle == WEST:
            self.x_minus = pos[0]+0.5
            self.y_minus = pos[1]
            self.x_plus = pos[0]-0.5
            self.y_plus = pos[1]

    def box_select(self, x0, y0, x1, y1):
        xs = [x0, x1]
        ys = [y0, y1]
        x,y = self.get_center_pos()

        if min(xs) <= x <= max(xs) and min(ys) <= y <= max(ys):
            self.force_select()

    def manual_place(self, event):
        self.init_create_component(event)

    def auto_place(self, auto_place_info):

        if self.x_minus == self.x_plus:
            # increasing y = SOUTH in tkinter
            if self.y_minus < self.y_plus:
                self.pos = [self.x_minus, (self.y_minus+self.y_plus)/2, SOUTH]
                self.create()
            else:
                self.pos = [self.x_minus, (self.y_minus+self.y_plus)/2, NORTH]
                self.create()
        elif self.y_minus == self.y_plus:
            if self.x_minus < self.x_plus:
                self.pos = [(self.x_minus+self.x_plus)/2, self.y_minus, EAST]
                self.create()
            else:
                self.pos = [(self.x_minus+self.x_plus)/2, self.y_minus, WEST]
                self.create()

    def request_value_label(self):
        window = RequestValueLabelWindow(self.canvas.master, self)
        self.canvas.master.wait_window(window)

    def import_tk_image(self):
        png = type(self).__name__
        if self.hover:
            png += '_hover'
        if self.selected:
            png += '_selected'
        png += '.png'

        if self.pos[2] is None:
            angle = self.init_angle
        else:
            angle = self.pos[2]

        img = Image.open(os.path.join(png_directory, png))
        size = round(self.canvas.grid_unit*(1-node_dot_radius))
        img = img.resize((size, size))
        img = img.rotate(angle)
        self.tk_image = ImageTk.PhotoImage(img)

    def create(self):
        self.add_or_replace_node_dots()
        x, y, angle = self.pos
        self.import_tk_image()
        self.image = self.canvas.create_image(
            *self.grid_to_canvas([x, y]), image=self.tk_image)
        self.add_or_replace_label()
        self.canvas.elements.append(self)
        self.set_allstate_bindings()

    def update_graphic(self):
        self.import_tk_image()
        self.canvas.itemconfig(self.image, image=self.tk_image)

    def adapt_to_grid_unit(self):
        self.update_graphic()
        self.canvas.coords(self.image, *self.grid_to_canvas(self.pos[:2]))
        self.add_or_replace_label()
        self.add_or_replace_node_dots()

    def double_click(self, event):
        self.modify_values(self)

    def modify_values(self, event=None):
        old_prop = self.prop
        self.request_value_label()
        if self.prop[0] is None and self.prop[1] is None:
            self.prop = old_prop
        else:
            self.add_or_replace_label()

    def init_create_component(self, event, angle=0.):
        if angle == 0. and  not self.canvas.used_arrows:
            self.canvas.message("Use arrows to rotate",t = 2)
        if angle != 0.:
            self.canvas.used_arrows = True

        self.canvas.track_changes = False
        self.init_angle = angle
        self.import_tk_image()
        self.image = self.canvas.create_image(
            self.canvas.canvasx(event.x), self.canvas.canvasy(event.y), image=self.tk_image)

        self.canvas.bind("<Button-1>", self.init_release)
        self.canvas.bind('<Motion>', self.init_on_motion)
        self.canvas.bind('<Escape>', self.abort_creation)
        self.canvas.bind(
            '<Left>', lambda event: self.init_create_component(event, angle=WEST))
        self.canvas.bind(
            '<Right>', lambda event: self.init_create_component(event, angle=EAST))
        self.canvas.bind(
            '<Up>', lambda event: self.init_create_component(event, angle=NORTH))
        self.canvas.bind(
            '<Down>', lambda event: self.init_create_component(event, angle=SOUTH))
        self.hover_enter(event)

    def unset_initialization_bindings(self):
        self.canvas.bind("<Button-1>", lambda event: None)
        self.canvas.bind('<Motion>', lambda event: None)
        self.canvas.bind('<Left>', lambda event: None)
        self.canvas.bind('<Right>', lambda event: None)
        self.canvas.bind('<Up>', lambda event: None)
        self.canvas.bind('<Down>', lambda event: None)
        self.canvas.bind('<Escape>', lambda event: None)

    def abort_creation(self, event=None):
        self.unset_initialization_bindings()
        self.canvas.delete(self.image)
        del self

    def init_release(self, event):
        self.unset_initialization_bindings()
        self.snap_to_grid()
        self.request_value_label()
        if self.prop[0] is None and self.prop[1] is None:
            self.abort_creation()
            return
        self.add_or_replace_label()
        self.canvas.elements.append(self)
        self.set_allstate_bindings()
        self.canvas.track_changes = True
        self.canvas.save()

    def open_right_click_menu(self, event):
        menu = tk.Menu(self.canvas, tearoff=0)
        menu.add_command(label="Edit", command=self.modify_values)
        menu.add_command(label="Rotate", command=self.rotate)
        menu.add_command(label="Delete", command=self.canvas.delete_selection)
        menu.add_separator()
        menu.add_command(label="Copy", command=self.canvas.copy_selection)
        menu.add_command(label="Cut", command=self.canvas.cut_selection)
        menu.tk_popup(event.x_root, event.y_root, 0)
        self.canvas.bind("<ButtonRelease-3>", lambda event: None)

    def rotate(self):
        if self.pos[2] % 180. == 0.:
            self.pos = [self.pos[0]-0.5, self.pos[1] +
                        0.5, (self.pos[2]+90.) % 360.]
        elif self.pos[2] % 180. == 90.:
            self.pos = [self.pos[0]+0.5, self.pos[1] -
                        0.5, (self.pos[2]+90.) % 360.]

        self.update_graphic()
        self.canvas.coords(self.image, *self.grid_to_canvas(self.pos[:2]))
        self.add_or_replace_label()
        self.add_or_replace_node_dots()
                
    def init_on_motion(self, event):
        x, y = self.canvas.coords(self.image)
        dx = self.canvas.canvasx(event.x) - x
        dy = self.canvas.canvasy(event.y) - y
        self.canvas.move(self.image, dx, dy)

    def move(self, dx, dy):
        '''
        Input given in canvas units
        '''
        self.canvas.move(self.image, dx, dy)
        if self.dot_minus is not None:
            self.canvas.delete(self.dot_minus)
            self.dot_minus = None
        if self.dot_plus is not None:
            self.canvas.delete(self.dot_plus)
            self.dot_plus = None
        self.add_or_replace_label()

    def on_updownleftright(self, event, angle):
        self.canvas.used_arrows = True
        gu = self.canvas.grid_unit
        if angle != self.pos[2]:
            self._angle = angle
            self.update_graphic()
            self.add_or_replace_label()

            if self.dot_minus is not None:
                self.canvas.delete(self.dot_minus)
                self.dot_minus = None
            if self.dot_plus is not None:
                self.canvas.delete(self.dot_plus)
                self.dot_plus = None
            self.was_rotated = True


    def delete(self, event=None):
        self.canvas.elements.remove(self)
        self.canvas.delete(self.image)
        if self.text is not None:
            self.canvas.delete(self.text)
        self.canvas.delete(self.dot_minus)
        self.canvas.delete(self.dot_plus)
        self.canvas.save()
        del self

    def add_or_replace_label(self):
        gu = self.canvas.grid_unit
        _, _, angle = self.pos
        x, y = self.canvas.coords(self.image)
        value, label = self.prop
        text = to_string(self.unit, label, value,
                         use_math=False, use_unicode=True)
        font = Font(family='Helvetica', size=int(gu/8.), weight='normal')
        text_position = (0.2)*gu
        if angle % 180. == 90. and self.text is None:
            self.text = self.canvas.create_text(
                x+text_position, y, text=text, anchor=tk.W, font=font)
        elif angle % 180. == 90. and self.text is not None:
            self.canvas.coords(self.text, x+text_position, y)
            self.canvas.itemconfig(self.text,
                                   text=text, anchor=tk.W, font=font)
        elif angle % 180. == 0. and self.text is None:
            self.text = self.canvas.create_text(
                x, y+text_position, text=text, anchor=tk.N, font=font)
        elif angle % 180. == 0. and self.text is not None:
            self.canvas.coords(self.text, x, y+text_position)
            self.canvas.itemconfig(self.text,
                                   text=text, anchor=tk.N, font=font)

    def snap_to_grid(self, event=None):
        x, y = self.canvas.coords(self.image)
        x0, y0 = self.canvas.canvas_center
        gu = float(self.canvas.grid_unit)

        if self.pos[2] is None:
            angle = self.init_angle
        else:
            angle = self.pos[2]

        if angle % 180. == 90.:
            self.pos = [
                round(float(x-x0)/gu),
                round(float(y-y0-gu/2.)/gu) + 0.5,
                angle]
            self.canvas.coords(self.image, *self.grid_to_canvas(self.pos[:2]))

        elif angle % 180. == 0.:
            self.pos = [
                0.5 + round(float(x-x0-gu/2.)/gu),
                round(float(y-y0)/gu),
                angle]
            self.canvas.coords(self.image, *self.grid_to_canvas(self.pos[:2]))

        # Add nodes in case of intersection with a wire
        self.add_nodes()

        # Add circles at the nodes of the component
        self.add_or_replace_node_dots()

class R(Component):
    """docstring for R"""

    def __init__(self, canvas, event=None, auto_place=None):
        self.unit = r'$\Omega$'
        super(R, self).__init__(canvas, event, auto_place)


class L(Component):
    """docstring for L"""

    def __init__(self, canvas, event=None, auto_place=None):
        self.unit = 'H'
        super(L, self).__init__(canvas, event, auto_place)


class C(Component):
    """docstring for C"""

    def __init__(self, canvas, event=None, auto_place=None):
        self.unit = 'F'
        super(C, self).__init__(canvas, event, auto_place)


class J(Component):
    """docstring for J"""

    def __init__(self, canvas, event=None, auto_place=None):
        self.unit = 'H'
        super(J, self).__init__(canvas, event, auto_place)


class G(Component):
    """docstring for J"""

    def __init__(self, canvas, event=None, auto_place=None):
        self.unit = ''
        super(G, self).__init__(canvas, event, auto_place)

    @property
    def prop(self):
        return [None, '']

    @prop.setter
    def prop(self, prop):
        pass

    def double_click(self, event):
        pass

    def request_value_label(self):
        pass

    def open_right_click_menu(self, event):
        menu = tk.Menu(self.canvas, tearoff=0)
        menu.add_command(label="Rotate", command=self.rotate)
        menu.add_command(label="Delete", command=self.delete)
        menu.add_separator()
        menu.add_command(label="Copy", command=self.canvas.copy_selection)
        menu.add_command(label="Cut", command=self.canvas.cut_selection)
        menu.tk_popup(event.x_root, event.y_root, 0)
        self.canvas.bind("<ButtonRelease-3>", lambda event: None)

    def add_or_replace_label(self):
        pass

    def add_or_replace_node_dots(self):
        gu = self.canvas.grid_unit
        canvas_coords_plus = self.grid_to_canvas([self.x_plus,self.y_plus])
        if self.dot_plus is None:
            self.dot_plus = self.canvas.create_circle(
                *canvas_coords_plus, gu*node_dot_radius)
        else:
            self.canvas.update_circle(self.dot_plus,
                *canvas_coords_plus, gu*node_dot_radius)

class RequestValueLabelWindow(tk.Toplevel):
    def __init__(self, master, component):
        tk.Toplevel.__init__(self, master)
        self.component = component

        # TODO add suggestions
        # TODO inform that filling two fields is optional
        fields = 'Value', 'Label'

        # Determine values of the fields
        v, l = self.component.prop
        if v is None:
            v = ''
        else:
            v = "%e" % v

        if l is None:
            l = ''

        field_values = [v, l]

        self.entries = []
        for i, field in enumerate(fields):
            row = tk.Frame(self)
            lab = tk.Label(row, width=7, text=field, anchor='w')
            ent = tk.Entry(row, width=7)
            ent.insert(tk.END, field_values[i])
            row.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)
            lab.pack(side=tk.LEFT)
            ent.pack(side=tk.RIGHT, expand=tk.YES, fill=tk.X)
            self.entries.append((field, ent))
        self.entries[0][1].focus()

        self.bind('<Return>', lambda event: self.ok())
        ok_button = tk.Button(self, text='OK', command=self.ok)
        ok_button.pack(side=tk.LEFT, padx=5, pady=5)
        cancel_button = tk.Button(self, text='Cancel', command=self.cancel)
        cancel_button.pack(side=tk.LEFT, padx=5, pady=5)

    def ok(self):
        value = self.entries[0][1].get()
        label = self.entries[1][1].get()
        if value.replace(' ', '') == "":
            v = None
        else:
            try:
                v = float(value)
            except ValueError:
                messagebox.showinfo(
                    "Incorrect value", "Enter a python style float, for example: 1e-2 or 0.01")
                self.focus_force()
                return None

        if label.replace(' ', '') == "":
            l = None
        else:
            l = label

        if l is None and v is None:
            messagebox.showinfo(
                "No inputs", "Enter a value or a label or both")
            self.focus_force()
            return None
        else:
            self.component.prop = [v, l]
            self.destroy()

    def cancel(self):
        self.destroy()


def open_canvas(netlist_file):
    # root = tk.Tk()
    # canvas = SnappingCanvas(root,
    #                         netlist_file=netlist_file, grid_unit=60, bg="white")
    # root.focus_force()
    # root.mainloop()
    app = MainWindow(tk.Tk(), netlist_file)
    app.mainloop()


class MainWindow(ttk.Frame):
    """ Main window class """

    def __init__(self, mainframe, netlist_file):
        """ Initialize the main Frame """
        ttk.Frame.__init__(self, master=mainframe)
        self.master.title('Circuit Editor')
        self.master.geometry('800x600')  # size of the main window
        self.master.rowconfigure(0, weight=1)  # make canvas expandable
        self.master.columnconfigure(0, weight=1)
        self.canvas = SnappingCanvas(
            self.master, netlist_file=netlist_file, grid_unit=60)


if __name__ == '__main__':
    open_canvas("test.txt")
