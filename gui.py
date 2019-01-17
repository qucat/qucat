try:
    # Tkinter for Python 2.xx
    import Tkinter as tk
    from tkFont import Font
except ImportError:
    # Tkinter for Python 3.xx
    import tkinter as tk
    from tkinter.font import Font
from PIL import Image, ImageTk
import numpy as np
import os
from bbq.utility import to_string
from copy import deepcopy

png_directory = os.path.join(os.path.dirname(__file__), ".graphics")


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


class SnappingCanvas(tk.Canvas):
    def __init__(self, master, grid_unit, netlist_file, **kw):
        tk.Canvas.__init__(self, master, bd=0, highlightthickness=0, **kw)
        self.netlist_file = netlist_file
        self.grid_unit = int(grid_unit)
        self.pack(fill=tk.BOTH, expand=1)
        self.focus_set()
        self.bind('r', lambda event: R(self, event))
        self.bind('l', lambda event: L(self, event))
        self.bind('c', lambda event: C(self, event))
        self.bind('j', lambda event: J(self, event))
        self.bind('w', lambda event: W(self, event))
        self.bind('s', lambda event: self.save())
        self.bind("<Configure>", self.draw_grid)
        self.bind('<Delete>', self.delete_selection)
        self.bind('<Control-a>', self.select_all)

        self.elements = []
        try:
            with open(netlist_file, 'r') as f:
                for el in f:
                    el = el.replace('\n', '')
                    el = el.split(";")
                    if el[0] in ['C', 'L', 'R', 'J', 'W']:
                        string_to_component(el[0], self, auto_place=el)

        except FileNotFoundError:
            with open(netlist_file, 'w') as f:
                pass

    def create_circle(self, x, y, r):
        x0 = x - r
        y0 = y - r
        x1 = x + r
        y1 = y + r
        return self.create_oval(x0, y0, x1, y1, fill='black')

    def draw_grid(self, event):
        self.delete("grid")
        dx = 1
        dy = 1
        w, h = event.width, event.height
        self.background = self.create_rectangle(
            0, 0, w, h, fill='white', tags='grid')
        for x in np.arange(self.grid_unit, w, self.grid_unit):
            for y in np.arange(self.grid_unit, h, self.grid_unit):
                self.create_line(x-dx, y, x+dx, y, tags='grid')
                self.create_line(x, y-dy, x, y+dy, tags='grid')
        self.tag_lower('grid')
        self.tag_bind('grid', '<ButtonPress-1>', self.start_selection_field)
        self.tag_bind('grid', "<B1-Motion>", self.expand_selection_field)
        self.tag_bind('grid', "<ButtonRelease-1>", self.end_selection_field)

    def start_selection_field(self, event):
        self.deselect_all()
        self.selection_rectangle = self.create_rectangle(
            event.x, event.y, event.x, event.y, dash=(3, 5))

    def expand_selection_field(self, event):
        self.deselect_all()
        x0, y0, x1, y1 = self.coords(self.selection_rectangle)
        self.coords(self.selection_rectangle, x0, y0, event.x, event.y)
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

    def delete_selection(self, event=None):
        for el in self.elements:
            if el.selected:
                el.delete()

    def save(self):
        # TODO auto save every x seconds
        with open(self.netlist_file, 'w') as f:
            for el in self.elements:
                if el.value is None:
                    v = ''
                else:
                    v = "%e" % el.value

                if el.label is None:
                    l = ''
                else:
                    l = el.label
                f.write("%s;%s;%s;%s;%s\n" % (
                    type(el).__name__,
                    el.coords_to_node_string(el.x_minus, el.y_minus),
                    el.coords_to_node_string(el.x_plus, el.y_plus),
                    v, l))


class TwoNodeElement(object):
    def __init__(self, canvas, event=None, auto_place=None):
        self.canvas = canvas
        self.grid_unit = canvas.grid_unit

        if auto_place is None and event is not None:
            self.x_center = event.x
            self.y_center = event.y
            self.manual_place(event)
        else:
            self.value = auto_place[3]
            self.label = auto_place[4]
            if self.label == '':
                self.label = None

            if self.value == '':
                self.value = None
            else:
                self.value = float(self.value)

            self.x_minus, self.y_minus = self.node_string_to_coords(
                auto_place[1])
            self.x_plus, self.y_plus = self.node_string_to_coords(
                auto_place[2])
            self.auto_place(auto_place)

    def coords_to_node_string(self, x, y):
        gu = self.grid_unit
        return "%d,%d" % (round(x/gu), round(y/gu))

    def node_string_to_coords(self, node):
        xy = node.split(',')
        x = int(xy[0])*self.grid_unit
        y = int(xy[1])*self.grid_unit
        return x, y

    def deselect(self):
        pass

    def force_select(self):
        pass

    def box_select(self, x0, y0, x1, y1):
        pass


class W(TwoNodeElement):
    def __init__(self, canvas, event=None, auto_place=None):
        self.value = None
        self.label = None
        self.hover = False
        self.selected = False
        super(W, self).__init__(canvas, event, auto_place)

    def manual_place(self, event):
        self.canvas.bind("<Button-1>", self.start_line)

    def auto_place(self, auto_place_info):
        self.create_component()

    def start_line(self, event):
        self.x_minus, self.y_minus = self.snap_to_grid(event)
        self.canvas.bind("<Motion>", self.show_line)
        self.canvas.bind("<Button-1>", self.end_line)

    def end_line(self, event):
        self.canvas.delete("temp")
        self.canvas.bind("<Button-1>", lambda event: None)
        self.canvas.bind('<Motion>', lambda event: None)

        self.x_plus, self.y_plus = self.snap_to_grid(event)
        self.coords_to_node_string(self.x_plus, self.y_plus)
        self.create_component()

    def create_component(self):
        self.line = self.canvas.create_line(
            self.x_minus, self.y_minus, self.x_plus, self.y_plus)
        self.dot_minus = self.canvas.create_circle(
            self.x_minus, self.y_minus, self.grid_unit/20.)
        self.dot_plus = self.canvas.create_circle(
            self.x_plus, self.y_plus, self.grid_unit/20.)
        self.canvas.elements.append(self)

    def show_line(self, event):
        self.canvas.delete("temp")
        self.canvas.create_line(
            self.x_minus, self.y_minus, event.x, event.y, tags='temp')

    def snap_to_grid(self, event):
        gu = float(self.grid_unit)
        return int(gu * round(float(event.x)/gu)),\
            int(gu * round(float(event.y)/gu))


class Component(TwoNodeElement):
    def __init__(self, canvas, event, auto_place):
        self.image = None
        self.value = None
        self.label = None
        self.hover = False
        self.selected = False
        self.text = None
        super(Component, self).__init__(canvas, event, auto_place)

    def manual_place(self, event):
        self.init_create_component(event)

    def auto_place(self, auto_place_info):

        if self.x_minus == self.x_plus:
            self.angle = -90.
            self.create_component(
                self.x_minus, (self.y_minus+self.y_plus)/2, self.angle)
        elif self.y_minus == self.y_plus:
            self.angle = 0
            self.create_component(
                (self.x_minus+self.x_plus)/2, self.y_minus, self.angle)
        self.add_label()
        self.canvas.elements.append(self)
        self.set_allstate_bindings()

    def request_value_label(self):
        window = RequestValueLabelWindow(self.canvas.master, self)
        self.canvas.master.wait_window(window)

    def import_tk_image(self, angle):
        png = type(self).__name__
        if self.hover:
            png += '_hover'
        if self.selected:
            png += '_selected'
        png += '.png'

        img = Image.open(os.path.join(png_directory, png))
        self.tk_image = ImageTk.PhotoImage(img.resize(
            (self.grid_unit, self.grid_unit)).rotate(angle))

    def create_component(self, x, y, angle=0.):
        self.x_center = x
        self.y_center = y
        self.angle = angle
        self.import_tk_image(angle)

        if self.image is not None:
            # Just replace tkimage
            self.canvas.itemconfig(self.image, image=self.tk_image)
        else:
            # Actually create image
            self.image = self.canvas.create_image(
                x, y, image=self.tk_image)

        if self.text is not None:
            self.add_label()

    def hover_enter(self, event):
        self.hover = True
        self.create_component(self.x_center, self.y_center, self.angle)
        self.canvas.tag_bind(self.image, "<Button-1>", self.on_click)
        self.canvas.tag_bind(self.image, "<B1-Motion>", self.on_motion)
        self.canvas.tag_bind(
            self.image, "<ButtonRelease-1>", self.release_motion)
        self.canvas.tag_bind(
            self.image, '<Double-Button-1>', self.modify_values)
        self.canvas.tag_bind(self.image, "<Shift-ButtonRelease-1>",
                             lambda event: self.release_motion(event, shift_control=True))
        self.canvas.tag_bind(self.image, "<Control-ButtonRelease-1>",
                             lambda event: self.release_motion(event, shift_control=True))

    def modify_values(self, event):
        self.request_value_label()
        self.add_label()

    def hover_leave(self, event):
        self.hover = False
        self.create_component(self.x_center, self.y_center, self.angle)

    def init_create_component(self, event, angle=0.):
        self.create_component(event.x, event.y, angle)
        self.canvas.bind("<Button-1>", self.init_release)
        self.canvas.bind('<Motion>', self.on_motion)
        self.canvas.bind(
            '<Left>', lambda event: self.init_create_component(event))
        self.canvas.bind(
            '<Right>', lambda event: self.init_create_component(event))
        self.canvas.bind(
            '<Up>', lambda event: self.init_create_component(event, angle=-90.))
        self.canvas.bind(
            '<Down>', lambda event: self.init_create_component(event, angle=-90.))

    def init_release(self, event):
        self.canvas.bind("<Button-1>", lambda event: None)
        self.canvas.bind('<Motion>', lambda event: None)
        self.canvas.bind('<Left>', lambda event: None)
        self.canvas.bind('<Right>', lambda event: None)
        self.canvas.bind('<Up>', lambda event: None)
        self.canvas.bind('<Down>', lambda event: None)

        self.snap_to_grid(event)
        self.request_value_label()
        self.add_label()
        self.canvas.elements.append(self)
        self.set_allstate_bindings()

    def set_allstate_bindings(self):
        self.canvas.tag_bind(self.image, "<Enter>", self.hover_enter)
        self.canvas.tag_bind(self.image, "<Leave>", self.hover_leave)

    def on_click(self, event):
        self.x_center_click = self.x_center
        self.y_center_click = self.y_center
        self.angle_click = self.angle
        self.x_center = event.x
        self.y_center = event.y

        self.canvas.bind('<Left>', lambda event: self.create_component(
            self.x_center, self.y_center))
        self.canvas.bind('<Right>', lambda event: self.create_component(
            self.x_center, self.y_center))
        self.canvas.bind('<Up>', lambda event: self.create_component(
            self.x_center, self.y_center, angle=-90.))
        self.canvas.bind('<Down>', lambda event: self.create_component(
            self.x_center, self.y_center, angle=-90.))

    def on_motion(self, event):
        dx = event.x - self.x_center
        dy = event.y - self.y_center
        self.canvas.move(self.image, dx, dy)
        if self.text is not None:
            self.canvas.move(self.text, dx, dy)
        self.x_center += dx
        self.y_center += dy

    def release_motion(self, event, shift_control=False):
        self.snap_to_grid(event)
        self.add_label()
        self.canvas.bind('<Left>', lambda event: None)
        self.canvas.bind('<Right>', lambda event: None)
        self.canvas.bind('<Up>', lambda event: None)
        self.canvas.bind('<Down>', lambda event: None)

        # if clicked without dragging or rotating
        if self.x_center_click == self.x_center \
                and self.y_center_click == self.y_center \
                and self.angle_click == self.angle:

            if shift_control:
                self.ctrl_shift_select()
            else:
                self.select()

    def select(self):
        self.canvas.deselect_all()
        if self.selected is False:
            self.selected = True
            self.create_component(self.x_center, self.y_center, self.angle)

    def box_select(self, x0, y0, x1, y1):
        xs = [x0, x1]
        ys = [y0, y1]
        if min(xs) <= self.x_center <= max(xs) and min(ys) <= self.y_center <= max(ys):
            self.force_select()

    def ctrl_shift_select(self):
        if self.selected is False:
            self.selected = True
            self.create_component(self.x_center, self.y_center, self.angle)
        elif self.selected is True:
            self.deselect()

    def force_select(self):
        self.selected = True
        self.create_component(self.x_center, self.y_center, self.angle)

    def deselect(self):
        self.selected = False
        self.create_component(self.x_center, self.y_center, self.angle)

    def delete(self, event=None):
        self.canvas.elements.remove(self)
        self.canvas.delete(self.image)
        if self.text is not None:
            self.canvas.delete(self.text)
        del self

    def add_label(self):

        if self.text is not None:
            self.canvas.delete(self.text)

        x, y = self.canvas.coords(self.image)
        text = to_string(self.unit, self.label, self.value,
                         use_math=False, use_unicode=True)
        font = Font(family='Helvetica', size=9, weight='normal')
        text_position = (0.3)*self.grid_unit
        if self.angle == -90.:
            self.text = self.canvas.create_text(
                x+text_position, y, text=text, anchor=tk.W, font=font)
        if self.angle == 0.:
            self.text = self.canvas.create_text(
                x, y+text_position, text=text, anchor=tk.N, font=font)

    def snap_to_grid(self, event):
        x, y = self.canvas.coords(self.image)
        gu = float(self.grid_unit)
        if self.angle == -90:
            self.x_center = int(gu * round(float(x)/gu))
            self.y_center = gu/2.+int(gu * round(float(y-gu/2.)/gu))
            self.canvas.coords(self.image, self.x_center, self.y_center)

            self.x_minus = self.x_center
            self.y_minus = self.y_center-gu/2.
            self.x_plus = self.x_center
            self.y_plus = self.y_center+gu/2.
        elif self.angle == 0.:
            self.x_center = gu/2.+int(gu * round(float(x-gu/2.)/gu))
            self.y_center = int(gu * round(float(y)/gu))
            self.canvas.coords(self.image, self.x_center, self.y_center)

            self.x_minus = self.x_center-gu/2.
            self.y_minus = self.y_center
            self.x_plus = self.x_center+gu/2.
            self.y_plus = self.y_center


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


class RequestValueLabelWindow(tk.Toplevel):
    def __init__(self, master, component):
        tk.Toplevel.__init__(self, master)
        self.component = component

        # TODO add suggestions
        # TODO inform that filling two fields is optional
        fields = 'Value', 'Label'

        # Determine values of the fields
        if component.value is None:
            v = ''
        else:
            v = "%e"%component.value
        if component.label is None:
            l = ''
        else:
            l = component.label
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
        if value == "":
            self.component.value = None
        else:
            self.component.value = float(value)

        if label == "":
            self.component.label = None
        else:
            self.component.label = str(label)
        self.destroy()

    def cancel(self):
        self.destroy()


def open_canvas(netlist_file):
    root = tk.Tk()
    canvas = SnappingCanvas(root,
                            netlist_file=netlist_file,
                            width=500, height=500, grid_unit=60, bg="white")
    root.focus_force()
    root.mainloop()


if __name__ == '__main__':
    open_canvas("test.txt")
