try:
    # Tkinter for Python 2.xx
    import Tkinter as tk
    from tkFont import Font
except ImportError:
    # Tkinter for Python 3.xx
    import tkinter as tk
    from tkinter.font import Font
from PIL import Image, ImageTk
import bbq.core_net
import numpy as np
bbq.core_net.pp={
    "element_width":1.,
    "element_height":1.,
    "margin":0.,
    "figsize_scaling":1.,
    "color":[0.15,0.15,0.15],
    "x_fig_margin":0.,
    "y_fig_margin":0.2,
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
        "lw":2
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

bbq.core_net.R(None,None,'').show(save_to = 'R.png',plot = False)
bbq.core_net.C(None,None,'').show(save_to = 'C.png',plot = False)
bbq.core_net.J(None,None,'').show(save_to = 'J.png',plot = False)
bbq.core_net.L(None,None,'').show(save_to = 'L.png',plot = False)

class SnappingCanvas(tk.Canvas):
    ''' A canvas that bites! ;-)'''
    def __init__(self, master, grid_unit, **kw):
        tk.Canvas.__init__(self, master,bd=0, highlightthickness=0, **kw)
        self.grid_unit = int(grid_unit)
        self.pack(fill=tk.BOTH, expand=1)
        self.focus_set()
        self.bind('r', lambda event: R(self,event))
        self.bind('l', lambda event: L(self,event))
        self.bind('c', lambda event: C(self,event))
        self.bind('j', lambda event: J(self,event))
        self.bind('w', lambda event: W(self,event))
        self.bind('s', lambda event: self.save())
        self.bind("<Configure>", self.draw_grid)
        self.elements = []

    def draw_grid(self,event):
        self.delete("grid")
        dx = 1
        dy = 1
        w, h = event.width, event.height
        for x in np.arange(self.grid_unit,w,self.grid_unit): 
            for y in np.arange(self.grid_unit,h,self.grid_unit):
                self.create_line(x-dx,y, x+dx,y, tags='grid')
                self.create_line(x,y-dy, x,y+dy, tags='grid')
        self.tag_lower('grid')

    def save(self):
        self.pretty_print_netlist(remove_wires = False)

    def pretty_print_netlist(self,remove_wires = False):
        for el in self.elements:
            # TODO: make sure %7s is enough space to write the node
            print("%-7s %-7s %s"%(el.node_minus,el.node_plus,el.comp.to_string(use_math = False)))
        print("")

class TwoNodeElement(object):
    def __init__(self, canvas, event):
        self.x = event.x
        self.y = event.y
        self.canvas = canvas
        self.grid_unit = canvas.grid_unit
        self.canvas.elements.append(self)

    def to_node(self,x,y):
        gu = self.grid_unit
        return "%d,%d"%(round(x/gu),round(y/gu))

    @property
    def node_minus(self):
        return self.comp.node_minus

    @node_minus.setter
    def node_minus(self,x,y):
        self.comp.node_minus = self.to_node(x,y)

    @property
    def node_plus(self):
        return self.comp.node_plus

    @node_plus.setter
    def node_plus(self,node):
        self.comp.node_plus = self.to_node(x,y)

class W(TwoNodeElement):
    def __init__(self, canvas, event):
        super(W, self).__init__(canvas, event)

        self.canvas.bind("<Button-1>", self.start_line)

    def start_line(self,event):
        self.x_start,self.y_start = self.snap_to_grid(event)
        self.canvas.bind("<Motion>", self.show_line)
        self.canvas.bind("<Button-1>", self.end_line)
        
    def end_line(self,event):
        self.canvas.delete("temp")
        self.canvas.bind("<Button-1>", lambda event: None)
        self.canvas.bind('<Motion>', lambda event: None)

        self.x_end,self.y_end = self.snap_to_grid(event)
        self.canvas.create_line(self.x_start, self.y_start,self.x_end,self.y_end)
        self.comp = bbq.core_net.W(self.to_node(self.x_start, self.y_start),
            self.to_node(self.x_end,self.y_end))

    def show_line(self,event):
        self.canvas.delete("temp")
        self.canvas.create_line(self.x_start, self.y_start,event.x,event.y,tags = 'temp')

    def snap_to_grid(self, event):
        gu = float(self.grid_unit)
        return int(gu * round(float(event.x)/gu)),\
                int(gu * round(float(event.y)/gu))
        
class Component(TwoNodeElement):
    def __init__(self, canvas, event):
        super(Component, self).__init__(canvas, event)
        self.image= None
        self.create_component(event)
    
    def request_value_label(self):
        window=RequestValueLabelWindow(self.canvas.master,self)
        self.canvas.master.wait_window(window)      

    def create_component(self,event,angle = 0.):
        self.angle = angle
        img = Image.open(self.png)
        self.tk_image = ImageTk.PhotoImage(img.resize(
            (self.grid_unit, self.grid_unit)).rotate(angle))
        if self.image is not None:
            self.canvas.delete(self.image)
        self.image= canvas.create_image(
            event.x,event.y, image=self.tk_image)

        self.canvas.bind("<Button-1>", self.init_release)
        self.canvas.bind('<Motion>', self.on_motion)
        self.canvas.bind('<Left>', lambda event: self.create_component(event))
        self.canvas.bind('<Right>', lambda event: self.create_component(event))
        self.canvas.bind('<Up>', lambda event: self.create_component(event, angle = -90.))
        self.canvas.bind('<Down>', lambda event: self.create_component(event, angle = -90.))
           
    def init_release(self,event):
        self.canvas.bind("<Button-1>", lambda event: None)
        self.canvas.bind('<Motion>', lambda event: None)
        self.canvas.bind('<Left>', lambda event: None)
        self.canvas.bind('<Right>', lambda event: None)
        self.canvas.bind('<Up>', lambda event: None)
        self.canvas.bind('<Down>', lambda event: None)

        self.value = None
        self.label = None
        x1,y1,x2,y2 = self.snap_to_grid(event)
        self.request_value_label()
        self.comp = self.component(
            self.to_node(x1,y1),
            self.to_node(x2,y2),
            self.value,self.label)
        self.add_label()

        self.canvas.tag_bind(self.image, "<Button-1>", self.on_click)
        self.canvas.tag_bind(self.image, "<B1-Motion>", self.on_motion)
        self.canvas.tag_bind(self.image, "<ButtonRelease-1>",self.snap_to_grid)

    

    def on_click(self, event):
        self.x = event.x
        self.y = event.y
        
    def on_motion(self, event):
        dx = event.x - self.x
        dy = event.y - self.y
        self.canvas.move(self.image,dx ,dy)
        self.x +=dx
        self.y +=dy

    def add_label(self):
        x,y = self.canvas.coords(self.image)
        text = self.comp.to_string(use_math = False,use_unicode = True)
        font = Font(family='Helvetica',size=9, weight='normal')
        text_position = (0.5-bbq.core_net.pp["y_fig_margin"])*self.grid_unit
        if self.angle == -90.:
            self.text = self.canvas.create_text(x+text_position,y,text = text,anchor = tk.W, font = font)
        if self.angle == 0.:
            self.text = self.canvas.create_text(x,y+text_position,text = text,anchor = tk.N, font = font)

    def snap_to_grid(self, event):
        x,y = self.canvas.coords(self.image)
        gu = float(self.grid_unit)
        if self.angle == -90:
            x_snap = int(gu * round(float(x)/gu))
            y_snap = gu/2.+int(gu * round(float(y-gu/2.)/gu))
            self.canvas.coords(self.image,x_snap,y_snap)
            return x_snap, y_snap-gu/2.,x_snap, y_snap+gu/2.
        elif self.angle == 0.:
            x_snap = gu/2.+int(gu * round(float(x-gu/2.)/gu))
            y_snap = int(gu * round(float(y)/gu))
            self.canvas.coords(self.image,x_snap,y_snap)
            return x_snap-gu/2., y_snap,x_snap+gu/2., y_snap


class R(Component):
    """docstring for R"""
    def __init__(self, canvas, event):
        self.png = 'R.png'
        self.component =bbq.core_net.R
        super(R, self).__init__(canvas, event)
class L(Component):
    """docstring for L"""
    def __init__(self, canvas, event):
        self.png = 'L.png'
        self.component =bbq.core_net.L
        super(L, self).__init__(canvas, event)
class C(Component):
    """docstring for C"""
    def __init__(self, canvas, event):
        self.png = 'C.png'
        self.component =bbq.core_net.C
        super(C, self).__init__(canvas, event)
class J(Component):
    """docstring for J"""
    def __init__(self, canvas, event):
        self.png = 'J.png'
        self.component =bbq.core_net.J
        super(J, self).__init__(canvas, event)

class RequestValueLabelWindow(tk.Toplevel):
    def __init__(self,master,component):
        tk.Toplevel.__init__(self, master)
        self.component = component

        # TODO add suggestions
        # TODO inform that filling two fields is optional
        fields = 'Value', 'Label'
        self.entries = []
        for field in fields:
          row = tk.Frame(self)
          lab = tk.Label(row, width=7, text=field, anchor='w')
          ent = tk.Entry(row, width=7)
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
            self.component.value  = None
        else:
            self.component.value  = float(value)

        if label == "":
            self.component.label  = None
        else:
            self.component.label  = str(label)
        self.destroy()

    def cancel(self):
        self.destroy()



        

root = tk.Tk()
canvas = SnappingCanvas(root, width=500, height=500,grid_unit = 60, bg="white")
root.focus_force()
root.mainloop()