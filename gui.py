try:
    # Tkinter for Python 2.xx
    import Tkinter as tk
except ImportError:
    # Tkinter for Python 3.xx
    import tkinter as tk
from bbq.core_sp  import R 
from PIL import Image, ImageTk
import bbq.core_sp
import numpy as np
bbq.core_sp.pp['element_width'] = 1.
bbq.core_sp.pp['element_height'] = 1.
bbq.core_sp.pp['margin'] = 0.
bbq.core_sp.pp['x_fig_margin'] = 0.
bbq.core_sp.pp['y_fig_margin'] = 0.2
bbq.core_sp.pp['R']['lw'] = 3
bbq.core_sp.pp['W']['lw'] = 3
R('').show(save_to = 'R.png',plot = False)

class SnappingCanvas(tk.Canvas):
    ''' A canvas that bites! ;-)'''
    def __init__(self, master, grid_unit, **kw):
        tk.Canvas.__init__(self, master,bd=0, highlightthickness=0, **kw)
        self.grid_unit = int(grid_unit)
        self.pack(fill=tk.BOTH, expand=1)
        self.focus_set()
        self.bind('r', lambda event: Component(self,event))
        self.bind('w', lambda event: Wire(self,event))
        self.bind('s', lambda event: save())
        self.bind("<Configure>", self.draw_grid)
        self.netlist = []

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
        print(self.netlist)

class Wire(object):
    def __init__(self, canvas, event):
        self.canvas = canvas
        self.grid_unit = canvas.grid_unit

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

    def show_line(self,event):
        self.canvas.delete("temp")
        self.canvas.create_line(self.x_start, self.y_start,event.x,event.y,tags = 'temp')

    def snap_to_grid(self, event):
        gu = float(self.grid_unit)
        return int(gu * round(float(event.x)/gu)),\
                int(gu * round(float(event.y)/gu))

class Component(object):
    def __init__(self, canvas, event):
        self.x = event.x
        self.y = event.y
        self.canvas = canvas
        self.grid_unit = canvas.grid_unit
        self.image= None
        self.create_component(event)
    
    def request_value_label(self):
        # TODO add suggestions
        # TODO inform that filling two fields is optional
        fields = 'Value', 'Label'

        def fetch(self,request_root,entries):
            self.value  = entries[0][1].get()
            self.label  = entries[1][1].get() 
            # TODO deal with no entry case
            request_root.quit()

        def makeform(request_root, fields):
           entries = []
           for field in fields:
              row = tk.Frame(request_root)
              lab = tk.Label(row, width=7, text=field, anchor='w')
              ent = tk.Entry(row, width=7)
              row.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)
              lab.pack(side=tk.LEFT)
              ent.pack(side=tk.RIGHT, expand=tk.YES, fill=tk.X)
              entries.append((field, ent))
           return entries

        def cancel(request_root):
            # TODO cancel creation of component
            request_root.quit()

        
        request_root = tk.Toplevel(root)
        entries = makeform(request_root, fields)
        request_root.bind('<Return>', (lambda event, e=entries: fetch(self,request_root,e)))   
        ok_button = tk.Button(request_root, text='OK', command=(lambda e=entries: fetch(self,request_root,e)))
        ok_button.pack(side=tk.LEFT, padx=5, pady=5)
        cancel_button = tk.Button(request_root, text='Cancel', command=cancel(request_root))
        cancel_button.pack(side=tk.LEFT, padx=5, pady=5)
        request_root.mainloop()

    def create_component(self,event,angle = 0.):
        self.angle = angle
        img = Image.open('R.png')
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

        self.request_value_label()
        self.on_release(event)

        self.canvas.tag_bind(self.image, "<Button-1>", self.on_click)
        self.canvas.tag_bind(self.image, "<B1-Motion>", self.on_motion)
        self.canvas.tag_bind(self.image, "<ButtonRelease-1>",self.on_release)

    def on_click(self, event):
        self.x = event.x
        self.y = event.y
        
    def on_motion(self, event):
        dx = event.x - self.x
        dy = event.y - self.y
        self.canvas.move(self.image,dx ,dy)
        self.x +=dx
        self.y +=dy

    def on_release(self, event):
        x,y = self.canvas.coords(self.image)
        gu = float(self.grid_unit)
        if self.angle == -90:
            self.canvas.coords(self.image,\
                int(gu * round(float(x)/gu)),\
                gu/2.+int(gu * round(float(y-gu/2.)/gu)))
        elif self.angle == 0.:
            self.canvas.coords(self.image,\
                gu/2.+int(gu * round(float(x-gu/2.)/gu)),\
                int(gu * round(float(y)/gu)))


    

        

root = tk.Tk()
canvas = SnappingCanvas(root, width=500, height=500,grid_unit = 60, bg="white")
root.focus_force()
root.mainloop()