try:
    # Tkinter for Python 2.xx
    import Tkinter as tk
    from tkFont import Font
    from Tkinter import tkMessageBox as messagebox
    from Tkinter import tkFileDialog as filedialog
except ImportError:
    # Tkinter for Python 3.xx
    import tkinter as tk
    from tkinter.font import Font
    from tkinter import messagebox,filedialog
from PIL import Image, ImageTk
from tkinter import ttk
import numpy as np
import os
from Qcircuits.src._utility import to_string
from Qcircuits.src._constants import *
from copy import deepcopy

png_directory = os.path.join(os.path.dirname(__file__), ".graphics")
node_dot_radius = 1./30.
lw = 1./50.
lw_hover = 2.*lw
lw_select_hover = 5.*lw
lw_select = 3.*lw


class CircuitEditor(tk.Canvas):
    def __init__(self, master, grid_unit, netlist_filename, **kw):
        """
        The CircuitEditor is the only widget which populates the MainWindow of the application.
        The main part of this widget is a Tkinter Canvas on which the user visualises
        and interacts with an electrical circuit.

        Circuit components are placed manually on the canvas and snap to the canvas grid.
        The user can navigate the canvas by scrolling to zoom in or out, or pan horizontally and vertically.
        The user can then drag and drop, edit, copy/cut/paste, etc... these components at will.
        Each change made by the user is automatically saved in a file defined by the netlist_filename.
        Each line of this text file is in the format:
        <type> (C,L,R...);<x,y (node_minus in grid unit)>;<x,y (node_plus in grid unit)>;value;symbol
        This file can be read by this application to load a circuit, or by an analysis software.
        These changes are also logged in a history variable which enables actions to be un/redone.

        Coordinate systems
        ==================

        The Canvas widget uses two coordinate systems; 
        the window coordinate system, with (0, 0) in the upper left corner
        and (canvas.winfo_width,canvas.winfo_height) in the lower right corner, 
        and a canvas coordinate system which specify where the items are drawn.

        By scrolling the canvas, you can specify which part of the canvas coordinate system 
        to show in the window.

        To convert from window coordinates to canvas coordinates, use the canvasx and canvasy methods.
        For example the upper left corner of the window has canvas coordinates (canvasx(0),canvasy(0))
        and the bottom right corner of the window has canvas coordinates 
        (canvasx(canvas.winfo_width), canvasy(canvas.winfo_height))

        A circuit component has thus a unique set of canvas coordinates (x_can, y_can).
        They can be seen on the window if:
        canvasx(0) < x_can < canvasx(canvas.winfo_width) and canvasx(0) < y_can < canvasy(canvas.winfo_height)

        We use an additional set of coordinates called grid coordinates, where the horizontal and
        vertical directions are divided into discrete steps to the circuit components will snap to.
        Any node of a circuit component can thus be assigned two (possible negative) integers, 
        (x_grid,y_grid) corresponding to a number of horizontal and vertical steps away from an 
        initially defined center. 
        The center and the size of the step size is defined in canvas units (canvas.canvas_center and
        canvas.grid_unit respectively) such that 
        x_can = canvas.canvas_center[0]+canvas.grid_unit*x_grid
        y_can = canvas.canvas_center[1]+canvas.grid_unit*y_grid

        canvas.canvas_center and canvas.grid_unit are not constant and are modified when zooming in/out

        Nodes, wires and intersections
        ==============================

        All circuit components currently implemented have two nodes except the special case
        of the ground element, which has one. 
        Nodes are represented as full black circles in the editor. 

        When a node (node 1) is placed on a wire, this wire automatically splits into two wires, each sharing
        a node with (node 1).
        Conversely, when a wire is moved, or is created such that it crosses a node (node 1), it will 
        split into two wires, each sharing a node with (node 1).
        
        However, if two wires cross, such that neither has a node placed on a wire, there will be no
        splitting of wires or creation of new nodes. This case implements a cros-over.
        """
        
        self.netlist_filename = netlist_filename
        '''In the netlist file is stored at all times 
        (except whilst drag/dropping elements)
        all the information necessary to re-construct 
        the circuit as the user sees it.'''

        self.grid_unit = int(grid_unit)
        '''Sets the size, in pixels of the 
        grid initially displayed upon opening 
        the editor. This variable'''

        self.elements = []
        '''List which stores all the circuit elements 
        currently placed on the canvas'''

        self.in_creation = None
        '''in_creation is None if no element is being created.
        It is equal to a circuit element if that 
        particular element is under creation.
        This allows us to handle zooming in/out 
        during the creation of a component'''

        self.copied_elements = []
        '''When elements are copied (or cut), 
        deepcopies of these elements are made and
        they are stored without being displayed 
        in this list'''

        self.history = []
        '''Everytime the user makes a change to the circuit, 
        this list is appended with a string representation of 
        all the components in the circuit
        hence allowing us to undo an arbitrary number of user actions '''

        self.history_location = -1
        '''Informs us of what index of the self.history variable
        the current circuit seen on the canvas corresponds'''

        self.track_changes = False
        '''We can cancel the appending of changes to the self.history
        variable (not to the netlist file)
        by setting this variable to True.
        This is useful for example when we want to delete all
        the components on the canvas: by calling undo afterwards
        we want to have them all reappear, not reappear one by one.
        Since deleting a single component is logged in the history 
        by default, we want to set track_changes to false, then
        delete all components, set track changes to true and finally
        save, which will add an entry to self.history.
        '''

        self.initialize_user_knowledge_tracking_variables()
        
        self.build_gridframe()
        self.build_menubar(master)
        self.build_scrollbars()
        self.build_canvas()
        
        self.configure_scrollbars()
        self.set_canvas_center()
        self.configure_canvas()
        
        self.set_keyboard_shortcuts_element_creation()
        self.set_keyboard_shortcuts_other()

        self.load_or_create_netlist_file()
        self.center_window_on_circuit()

        # start tracking changes and save once
        # to initialize self.history
        self.track_changes = True
        self.save()

        # Sets the "active window" in your OS, 
        # the one towards which key strokes will be 
        # directed to be this circuit editor
        self.focus_set()
        
    ###########################
    # Initialization functions
    ###########################
    
    def build_menubar(self,master):
        '''
        Builds the File, Edit, ... menu bar situated at the top of
        the window.
        '''

        # initialize the menubar object
        self.menubar = tk.Menu(self.frame)

        ####################################
        # Define the label formatting
        ####################################

        # File, Edit, ... are defined to have a width of 6 characters
        menu_label_template = "{:<6}"

        # The items appearing in the cascade menu that appears when 
        # clicking on File for example will have 15 characters width on the 
        # left where the name of the functionality is provided and
        # 6 characters on the right where the keyboard shortcut is provided
        label_template = "{:<15}{:>6}"

        # The font for the cascading items is defined
        # Note: for everything to be aligned the chosen font should 
        # have a constant character width, the only one satisfying this 
        # condition is Courier new.
        #TODO use images as cascade menu items with an aligned and pretty font
        menu_font = Font(family="Courier New", size=9, weight='normal')

        ####################################
        # FILE cascade menu build
        ####################################

        # add new item to the menubar
        menu = tk.Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(
            label=menu_label_template.format("File"), 
            menu=menu)

        # add cascade menu items
        menu.add_command(
            label=label_template.format("Open", ""), 
            command=self.file_open, 
            font=menu_font)
        menu.add_command(
            label=label_template.format("Save", "Ctrl+S"), 
            command=self.save, 
            font=menu_font)
        menu.add_command(
            label=label_template.format("Exit", ""), 
            command=master.destroy, 
            font=menu_font)

        ####################################
        # EDIT cascade menu build
        ####################################

        # add new item to the menubar
        menu = tk.Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(
            label=menu_label_template.format("Edit"), 
            menu=menu)

        # add cascade menu items
        menu.add_command(
            label=label_template.format("Undo", "Ctrl+Z"), 
            command=self.ctrl_z, 
            font=menu_font)
        menu.add_command(
            label=label_template.format("Redo", "Ctrl+Y"), 
            command=self.ctrl_y, 
            font=menu_font)
        menu.add_separator()
        menu.add_command(
            label=label_template.format("Cut", "Ctrl+X"), 
            command=self.cut_selection, 
            font=menu_font)
        menu.add_command(
            label=label_template.format("Copy", "Ctrl+C"), 
            command=self.copy_selection, 
            font=menu_font)
        menu.add_command(
            label=label_template.format("Paste", "Ctrl+V"), 
            command=(lambda :self.event_generate('<Control-v>')), 
            font=menu_font)
        menu.add_separator()
        menu.add_command(
            label=label_template.format("Select all", "Ctrl+A"), 
            command=self.select_all, 
            font=menu_font)
        menu.add_separator()
        menu.add_command(
            label=label_template.format("Delete", "Del"), 
            command=self.delete_selection, 
            font=menu_font)
        menu.add_command(
            label=label_template.format("Delete all", ""), 
            command=self.delete_all, 
            font=menu_font)


        ####################################
        # INSERT cascade menu build
        ####################################
        
        # add new item to the menubar
        menu = tk.Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(
            label=menu_label_template.format("Insert"), 
            menu=menu)

        # add cascade menu items
        menu.add_command(
            label=label_template.format("Wire", "W"), 
            command=(lambda: self.event_generate('w')), 
            font=menu_font)
        menu.add_command(
            label=label_template.format("Junction", "J"), 
            command=(lambda: self.event_generate('j')), 
            font=menu_font)
        menu.add_command(
            label=label_template.format("Inductor", "L"), 
            command=(lambda: self.event_generate('l')), 
            font=menu_font)
        menu.add_command(
            label=label_template.format("Capacitor", "C"), 
            command=(lambda: self.event_generate('c')), 
            font=menu_font)
        menu.add_command(
            label=label_template.format("Resistor", "R"), 
            command=(lambda: self.event_generate('r')), 
            font=menu_font)
        menu.add_command(
            label=label_template.format("Ground", "G"), 
            command=(lambda: self.event_generate('g')), 
            font=menu_font)

        
        ####################################
        # VIEW cascade menu build
        ####################################
        
        # add new item to the menubar
        menu = tk.Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(
            label=menu_label_template.format("View"), 
            menu=menu)

        # add cascade menu items
        menu.add_command(
            label=label_template.format("Zoom in", "Ctrl+scroll"), 
            command=(lambda: self.zoom('in')), 
            font=menu_font)
        menu.add_command(
            label=label_template.format("Zoom out", "Ctrl+scroll"), 
            command=(lambda: self.zoom('out')), 
            font=menu_font)
        menu.add_separator()
        menu.add_command(
            label=label_template.format("Re-center", ""), 
            command=(self.center_window_on_circuit), 
            font=menu_font)

        # Add the menubar to the application
        master.config(menu=self.menubar)
    def build_gridframe(self):
        '''
        Builds the main area of the window (called a frame), 
        which should stick to the edges of the window and 
        expand as a user expands the window.

        This frame will be divided into a grid hosting the
        canvas, menubar, scrollbars
        '''

        # Builds a new frame, which will be divided into a grid
        # hosting the canvas, menubar, scrollbars
        self.frame = ttk.Frame()

        # Places the Frame widget self.frame in the parent 
        # widget (MainWindow) in a grid
        self.frame.grid()  

        # Configure the frames grid
        self.frame.grid(sticky='nswe')  # make frame container sticky
        self.frame.rowconfigure(0, weight=1)  # make canvas expandable in x
        self.frame.columnconfigure(0, weight=1)  # make canvas expandable in y
    def build_scrollbars(self):
        '''
        Builds horizontal and vertical scrollbars and places
        them in the window
        '''
        
        # Vertical and horizontal scrollbars for canvas
        self.hbar = ttk.Scrollbar(self.frame, orient='horizontal')
        self.vbar = ttk.Scrollbar(self.frame, orient='vertical')
        self.hbar.grid(row=1, column=0, sticky='we')
        self.vbar.grid(row=0, column=1, sticky='ns')
    def build_canvas(self):
        '''
        Initializes the canvas from which this object inherits and 
        places it in the grid of our window

        Actually in the previously called build_* functions, we have
        defined frames, menubars, etc.. which are seperate from the canvas.
        It is however convinient to do so since most of these definied buttons
        or scrollbars will be acting on the canvas itself.
        '''
        tk.Canvas.__init__(
            self, 
            self.frame, 
            bd=0, 
            highlightthickness=0,
            xscrollcommand=self.hbar.set, 
            yscrollcommand=self.vbar.set, 
            confine=False, 
            bg="white")
        self.grid(row=0, column=0, sticky='nswe')
    def configure_scrollbars(self):
        '''
        Define what functions the scrollbars should call
        when we interact with them.
        '''
        self.hbar.configure(command=self.scroll_x)  # bind scrollbars to the canvas
        self.vbar.configure(command=self.scroll_y)
    def set_canvas_center(self):
        '''
        Calculate the center of the canvas in 
        canvas coordinates. This information is necessary to 
        convert grid units to canvas units

        '''
        # Wait for the canvas to pop up before asking 
        # for its window size below
        self.update()  

        # winfo_width/height gives the 
        # width/height of the window in window units
        # canvasx/y converts window units to 
        # canvas units
        self.canvas_center = [
            self.canvasx(self.winfo_width()/2.),
            self.canvasy(self.winfo_height()/2.)]
    def configure_canvas(self):
        '''
        Ensure that the function "on_resize"
        is called each time the user resizes
        the window. 
        '''
        self.bind("<Configure>", self.on_resize)
    def set_keyboard_shortcuts_element_creation(self):
        '''
            Assign key R to the creation of a resistor, 
            C to a capacitor, etc...
        '''
        self.bind('r', lambda event: R(self, event))
        self.bind('l', lambda event: L(self, event))
        self.bind('c', lambda event: C(self, event))
        self.bind('j', lambda event: J(self, event))
        self.bind('w', lambda event: W(self, event))
        self.bind('g', lambda event: G(self, event))
    def set_keyboard_shortcuts_other(self):
        '''
        Assign keystrokes to functionalities
        accessible in the FILE and EDIT menus, 
        as well as configure what happens when the 
        user scrolls in combination with CTRL/SHIFT.
        '''
        
        #############################
        # FILE menu functionalities
        #############################
        self.bind('<Control-s>', self.save)

        #############################
        # EDIT menu functionalities
        #############################
        self.bind('<Delete>', self.delete_selection)
        self.bind('<Control-c>', self.copy_selection)
        self.bind('<Control-x>', self.cut_selection)
        self.bind('<Control-v>', self.paste)
        self.bind('<Control-a>', self.select_all)
        self.bind('<Control-y>', self.ctrl_y)
        self.bind('<Control-z>', self.ctrl_z)

        

        #############################
        # Mouse wheel functionalities
        #############################

        # zoom for Windows and MacOS, but not Linux
        self.bind('<Control-MouseWheel>', self.scroll_zoom)
        # zoom for Linux, wheel scroll down
        self.bind('<Control-Button-5>',   self.scroll_zoom)
        # zoom for Linux, wheel scroll up
        self.bind('<Control-Button-4>',   self.scroll_zoom)
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
    def load_or_create_netlist_file(self):
        '''
        If the file called "self.netlist_filename" is found, load it, if not, create 
        a blank one and load that.
        '''

        try:
            with open(self.netlist_filename, 'r') as f:
                netlist_file_string = [line for line in f]
        except FileNotFoundError:
            netlist_file_string = []
            with open(self.netlist_filename, 'w') as f:
                pass
        self.load_netlist(netlist_file_string)
    def initialize_user_knowledge_tracking_variables(self):
        '''
        Initialize a set of variables which will allow us
        track what functionalities the user has employed so far.
        This will allow us to provide hints telling the user 
        how to do stuff that he hasnt done yet.
        '''
        
        self.used_arrows = False
        '''If the user uses arrows to rotate an 
        element whilst creating it, we set this variable 
        to True, and hence stop hinting that he can do that.
        '''

    #############################
    #  CORE functions
    ##############################

    def elements_list_to_netlist_string(self):
        '''
        Maps the list of elements (self.elements)
        to a string which can be saved in 
        a file and read either by the gui code
        or by the Qcircuits code.

        Returns
        --------
        netlist_string: string
            of the form 
            "C;0,1;0,0;1.000000e-13;
            J;1,1;1,0;;L_J
            etc..."

        '''
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
        return netlist_string

    def save(self, event=None):
        '''
            Save the current state of the circuit to the file.
            If we are tracking changes (self.track_changes == True)
            also append the history list.
        '''

        # Obtain the string representation 
        # of all the elements on the canvas       
        netlist_string = self.elements_list_to_netlist_string()

        # Write that netlist string to the file
        with open(self.netlist_filename, 'w') as f:
            f.write(netlist_string)

        # If we are tracking the changes made to the circuit
        if self.track_changes:
            
            # In case we just used undo (ctrl-z), 
            # we first want to delete the all entries 
            # in history which are after the current state of the circuit
            del self.history[self.history_location+1:]

            # we append all the information about the current circuit
            # to the history list
            self.history.append(netlist_string)

            # and increase our location in the history 
            # lsit by one
            self.history_location += 1

        # Inform the user that the circuit was just saved
        self.message("Saving...")
    
    def load_netlist(self, lines):
        '''
        Loads a circuit stored as a list of strings, each string
        representing a component of the circuit. The format of these
        strings is the same as that of the lines in the netlist file

        Parameters
        ----------
        node_minus: list of strings
                    Each string should be in the format
                    <type> (C,L,R...);<x,y (node_minus in grid unit)>;<x,y (node_plus in grid unit)>;value;symbol
        '''

        self.delete_all(track_changes=False)
        for el in lines:
            el = el.replace('\n', '')
            el = el.split(";")

            #TODO replace this with 
            # less verbose code?
            if el[0] == 'W':
                W(self, auto_place=el)
            elif el[0] == 'R':
                R(self, auto_place=el)
            elif el[0] == 'L':
                L(self, auto_place=el)
            elif el[0] == 'J':
                J(self, auto_place=el)
            elif el[0] == 'C':
                C(self, auto_place=el)
            elif el[0] == 'G':
                G(self, auto_place=el)

    def on_resize(self, event=None):
        '''
        Called when the user resizes the window, see configure_scrollregion
        and draw_grid for more detail
        '''
        self.configure_scrollregion()
        self.draw_grid(event)

    def draw_grid(self, event=None):
        '''
        Called when the user resizes the window. 
        Will delete and rebuild the grid.
        '''

        # Delete old grid
        self.delete("grid")

        # Get visible area of the canvas in canvas units
        box_canvas = (self.canvasx(0),  
                      self.canvasy(0),
                      self.canvasx(self.winfo_width()),
                      self.canvasy(self.winfo_height()))

        # Create a white background
        # it is invisible, but we can assign actions to the user clicking on it
        self.background = self.create_rectangle(
            *box_canvas, fill='white', outline='', tags='grid')

        # Create x and y coordinates for the grid
        grid_x = np.arange(
            self.canvas_center[0], box_canvas[2], self.grid_unit).tolist()
        grid_x += np.arange(self.canvas_center[0]-self.grid_unit,
                            box_canvas[0], -self.grid_unit).tolist()
        grid_y = np.arange(
            self.canvas_center[1], box_canvas[3], self.grid_unit).tolist()
        grid_y += np.arange(self.canvas_center[1]-self.grid_unit,
                            box_canvas[1], -self.grid_unit).tolist()

        # write the grid lines
        for x in grid_x:
            for y in grid_y:
                self.create_line(x-1, y, x+2, y, tags='grid')
                self.create_line(x, y-1, x, y+2, tags='grid')

        # Put the grid behind all other elements of the canvas
        self.tag_lower('grid')

        # Determine what happens when the user clicks on the background
        self.tag_bind('grid', '<ButtonPress-1>', self.start_selection_field)
        self.tag_bind('grid', "<B1-Motion>", self.expand_selection_field)
        self.tag_bind('grid', "<ButtonRelease-1>", self.end_selection_field)
        self.tag_bind('grid', "<Button-3>", self.right_click)

    ###########################
    # FILE menu functionalities
    ###########################

    def file_open(self):
        '''
        Triggered by the menu bar button File>Open.
        Opens a dialog window where the user can choose a file
        then loads the file.
        '''

        # Prompt user for file name
        netlist_filename = filedialog.askopenfilename(initialdir = os.getcwd())

        if netlist_filename == '':
            # User cancelled
            pass
        else:

            # open file
            with open(netlist_filename, 'r') as f:

                # extract the netlist file string
                netlist_file_string = [line for line in f]
            
            # here we don't want to track changes
            # as we're creating many components in one go, 
            # we'll be saving the new circuit afterwards
            self.track_changes = False
            try:
                # try and load the netlist
                self.load_netlist(netlist_file_string)
            except Exception as e:
                # in case the content of the file was not in the right format
                self.message("Not all components of the file could be loaded")
                print("Loading file failed with error:")
                print(e)
            else:
                self.message("File succesfully loaded")
            
            # Save changes and turn change tracker back on
            self.track_changes = True
            self.save()
            
            # Center window in case the other circuit was 
            # built at a different location on the canvas
            self.center_window_on_circuit()

    #############################
    # SCROLLING/ZOOMING
    ##############################

    def scroll_y_wheel(self, event):
        '''
        Triggered by the user scrolling (in combination with no particular key presses).
        '''

        # Determine which direction the user is scrolling 
        # if using windows, then event.delta has also a different
        # amplitude depending on how fast the user is scrolling, 
        # but we ignore that 
        if event.num == 5 or event.delta < 0:
            direction = 1
        if event.num == 4 or event.delta > 0:
            direction = -1

        # Move the canvas appropriately
        self.yview_scroll(direction, tk.UNITS)

        # reconfigure the region in which the scroll bars
        # can scroll
        self.configure_scrollregion()

        # redraw the grid such that it fills the 
        # visible canvas
        self.draw_grid(event)

    def scroll_x_wheel(self, event):
        '''
        Triggered by the user is SHIFT+scrolling 
        '''

        # Determine which direction the user is scrolling 
        # if using windows, then event.delta has also a different
        # amplitude depending on how fast the user is scrolling, 
        # but we ignore that.
        # Note: Linux -> event.num and  Windows -> event.delta
        if event.num == 5 or event.delta < 0:
            direction = 1
        if event.num == 4 or event.delta > 0:
            direction = -1

        # Move the canvas appropriately
        self.xview_scroll(direction, tk.UNITS)

        # reconfigure the region in which the scroll bars
        # can scroll
        self.configure_scrollregion()

        # redraw the grid such that it fills the 
        # visible canvas
        self.draw_grid(event)
    
    def scroll_zoom(self, event):
        '''
        Called when the user ALT+Scrolls.
        Zooms in/out of the canvas.
        Zooming works by chaning the grid_unit of the canvas
        and re-plotting all the circuit elements.
        During zooming, we move the circuit and grid such that
        the position of the mouse on the grid remains constant
        '''
        
        # Sets the smallest/largest allowed grid_unit size
        smallest_grid_unit = 35
        largest_grid_unit = 100

        # Determine which direction the user is scrolling 
        # if using windows, then event.delta has also a different
        # amplitude depending on how fast the user is scrolling, 
        # but we ignore that.
        # Note: Linux -> event.num and  Windows -> event.delta
        old_grid_unit = self.grid_unit
        
        # Default scaling of the grid_unit for slow scrolling on windows
        # or all scrolling on other OS
        scaling = 1.08

        # If on windows, we can change the scaling in case of fast scrolling
        try:
            if abs(event.delta) > 120:
                scaling = 1.15
        except:
            pass

        # Determine which direction the user is scrolling 
        # and scale the grid_unit accordingly
        # Note: Linux -> event.num and  Windows -> event.delta
        if event.num == 5 or event.delta < 0:  # scroll out, smaller
            new_grid_unit = int(self.grid_unit/scaling)
            if new_grid_unit == old_grid_unit:
                new_grid_unit -= 1
        elif event.num == 4 or event.delta > 0:  # scroll in, bigger
            new_grid_unit = int(self.grid_unit*scaling)
            if new_grid_unit == old_grid_unit:
                new_grid_unit += 1

        # If the user is trying to go below/above the 
        # smallest/largest grid unit size, just
        # display a message
        if smallest_grid_unit > new_grid_unit:
            self.message("Can't zoom out more")
        elif new_grid_unit > largest_grid_unit:
            self.message("Can't zoom in more")

        else:

            # position of the mouse when the scrolling occured
            # in old grid units
            grid_mouse_pos_old = self.canvas_to_grid(
                [self.canvasx(event.x), self.canvasy(event.y)])

            # change the grid unit
            self.grid_unit = new_grid_unit

            
            # position of the mouse when the scrolling occured
            # in canvas units
            canvas_mouse_pos = self.grid_to_canvas(grid_mouse_pos_old)

            # Amount we have to shift the canvas such that the 
            # position of the mouse on the new and old grid remains constant
            canvas_center_shift = [
                self.canvasx(event.x)-canvas_mouse_pos[0], 
                self.canvasy(event.y)-canvas_mouse_pos[1]]

            # Shift the center of the canvas such that 
            # position of the mouse on the new and old grid remains constant
            self.canvas_center = [self.canvas_center[0]+canvas_center_shift[0],
                                  self.canvas_center[1]+canvas_center_shift[1]]

            # Move and scale all ALREADY CREATED elements to adapt to the 
            # new grid
            for el in self.elements:
                el.adapt_to_grid_unit()

            # Move and scale all IN CREATION elements to adapt to the 
            # new grid    
            if self.in_creation is not None:
                self.in_creation.init_adapt_to_grid_unit(event)

            # redraw the grid such that it fills the 
            # visible canvas
            self.draw_grid(event)

            # reconfigure the region in which the scroll bars
            # can scroll
            self.configure_scrollregion()

       
    def configure_scrollregion(self):

        '''
        Called every time some moving around or zooming occurs on the canvas
        Configures the range that is scrollable with the scrollbars
        '''

        extra_scrollable_region = 50 # in canvas units
        '''
        if the circuit fulls the visible canvas, there will 
        still be a small gap in the scrollbar to indicate to
        the user he can use the scrollbars to scroll down by some 
        small amount
        '''

        # get visible area of the canvas in canvas units
        box_canvas = [self.canvasx(0)-extra_scrollable_region,  
                      self.canvasy(0)-extra_scrollable_region,
                      self.canvasx(self.winfo_width())+extra_scrollable_region,
                      self.canvasy(self.winfo_height())+extra_scrollable_region]

        # If there are some drawn circuit elemnts
        # set box_elemnts to describe the area filled by the circuit 
        # in canvas units
        if len(self.elements) > 0:
            xs = [el.x_minus for el in self.elements] + \
                [el.x_plus for el in self.elements]
            ys = [el.y_minus for el in self.elements] + \
                [el.y_plus for el in self.elements]
            box_elements = self.grid_to_canvas(
                [min(xs)-1, min(ys)-1])+self.grid_to_canvas([max(xs)+1, max(ys)+1])

            # If there are some drawn circuit elemnts, the scrollable region
            # should show that the user has some circuit elements to discover
            # if he scrolls down/up/left/right
            self.configure(
                scrollregion=[min(box_elements[0], box_canvas[0]), min(box_elements[1], box_canvas[1]),
                              max(box_elements[2], box_canvas[2]), max(box_elements[3], box_canvas[3])])
        else:
            # if there are no drawn circuit elements, just indicate
            # that the user can scroll down/up/left/right a little bit
            self.configure(scrollregion=box_canvas)

    def zoom(self, direction = 'in'):
        '''
        Generates an event which simulates the user
        scrolling by one increment

        Parameters
        ----------
        direction:  string
                    'in' tp scroll in
                    'out' to scroll out
        '''

        # Location at which the fake scrolling occurs
        kwargs = {'x':0,'y':0}

        args = ['<Control-MouseWheel>']
        if direction == 'in':
            self.event_generate(*args, delta = 121, **kwargs)
        if direction == 'out':
            self.event_generate(*args, delta = -121, **kwargs)

    def scroll_x(self, *args, **kwargs):
        """ 
        Is called when the user interacts with the horizontal scroll bar
        """
        # shift canvas horizontally
        self.xview(*args)

        # redraw the grid such that it fills the 
        # visible canvas
        self.draw_grid()

    def scroll_y(self, *args, **kwargs):
        """ 
        Is called when the user interacts with the vertical scroll bar
        """
        # shift canvas vertically
        self.yview(*args)

        # redraw the grid such that it fills the 
        # visible canvas
        self.draw_grid()

    #############################
    #  COPY/CUT/PASTE
    ##############################

    def cut_selection(self, event=None):
        '''
        Called on CTRL+X or Edit>Cut.
        Copies and deletes all selected elements.
        '''
        self.copy_selection()
        self.delete_selection()

    def copy_selection(self, event=None):
        '''
        Called on CTRL+C or Edit>Copy.
        Adds selected elements to the self.copied_elements variable.
        '''
        
        # Since deepcopying calls the __init__ of elements
        # we forbid any additions the history variable to be on the safe side
        self.track_changes = False

        self.copied_elements = [deepcopy(el)
                                for el in self.elements if el.selected]
        self.track_changes = True

    def paste(self, event=None):
        '''
        Called on CTRL+V or Edit>Paste.
        Creates all the elements in self.copied_elements
        and has them hover under the users mouse until
        he clicks somewhere, snapping them in place on the grid.
        '''

        if len(self.copied_elements) > 0:
            self.deselect_all()

            # Create a list of elements to paste, this allows
            # us to paste multiple times the same content
            # without the different version of the elements having
            # anything in common.
            self.track_changes = False
            to_paste = [deepcopy(el) for el in self.copied_elements]
            self.track_changes = True

            # Calculate top left position of the circuit to paste
            # in grid units then convert it to canvas units
            x_min = min([el.x_minus for el in to_paste] +
                        [el.x_plus for el in to_paste])
            y_min = min([el.y_minus for el in to_paste] +
                        [el.y_plus for el in to_paste])
            x_min, y_min = self.grid_to_canvas([x_min, y_min])

            # shift in position to apply to all elements
            # such that the top left of the circuit lies
            # under the mouse
            # in canvas units
            dx = self.canvasx(event.x)-x_min
            dy = self.canvasy(event.y)-y_min

            for el in to_paste:

                # create the component
                el.create()
                # scale the component and position it where it was copied
                el.adapt_to_grid_unit()
                # select it
                el.force_select()
                # move it such that the top left of the circuit lies
                # under the mouse
                el.move(dx, dy)

                # used to update the position of components label
                el.add_or_replace_label()

            #########################
            # Note: all the bindings below need only
            # be applied to a single element since
            # "on_motion", "release_motion_paste" acts
            # on all the selected elements
            # and we have set all the pasted elements to 
            # be selected above
            #########################

            # Ensure that when the mouse moves, the selection 
            # moves such that the top left of the circuit lies
            # under the mouse
            self.bind("<Motion>", el.on_motion)
            
            if len(to_paste) == 1:
                
                # Snaps the element to the grid and removes the binding 
                # of arrow keys and replaces binding of button-press to box selection
                self.bind("<ButtonPress-1>", el.release_motion_paste_single)

                # If only a single element is pasted, allow the user to rotate that 
                # element with arrows
                self.bind('<Left>', lambda event: el.on_updownleftright(event, angle=WEST))
                self.bind('<Right>', lambda event: el.on_updownleftright(event, angle=EAST))
                self.bind('<Up>', lambda event: el.on_updownleftright(event, angle=NORTH))
                self.bind('<Down>', lambda event: el.on_updownleftright(event, angle=SOUTH))
            else:
                # Snaps the elements (the selection) to the grid and
                # replaces binding of button-press to box selection
                self.bind("<ButtonPress-1>", el.release_motion_paste)


    #############################
    #  HISTORY MANAGEMENT
    ##############################

    def ctrl_z(self, event=None):
        '''
        Called on CTRL+Z or Edit>Undo.
        Returns to the circuit configuration previously stored in the
        self.history variable. 
        A circuit is added to the self.history variable when self.save()
        is called and self.track_changes == True.
        '''

        if self.history_location > 0:
            self.track_changes = False
            self.history_location -= 1
            self.load_netlist(self.history[self.history_location].split('\n'))
            self.save()
            self.track_changes = True
        else:
            self.message('Nothing to undo')

    def ctrl_y(self, event=None):
        '''
        Called on CTRL+Y or Edit>Redo.
        Returns to the next circuit configuration available in the
        self.history variable. 
        This is possible if the user just undid (CTRL-Z) and did 
        not yet create or move a component.
        Indeed if the user does the latter then self.save() should be
        called with self.track_changes == True, which erases all future 
        self.history entries.
        This function will then display the message 'Nothing to redo'.
        '''

        if 0 <= self.history_location < len(self.history)-1:
            self.track_changes = False
            self.history_location += 1
            self.load_netlist(self.history[self.history_location].split('\n'))
            self.track_changes = True
        else:
            self.message('Nothing to redo')

    #############################
    #  RIGHT CLICK
    ##############################

    def right_click(self, event):
        '''
        Called when the user right clicks on the grid or on the 
        background of the canvas.
        Deselects all components and opens the right click menu.
        '''
        self.deselect_all()
        self.bind("<ButtonRelease-3>", self.open_right_click_menu)

    def open_right_click_menu(self, event):
        '''
        Menu opened when the user right clicks on the grid or on the 
        background of the canvas.
        Deselects all components and opens the right click menu.
        '''
        menu = tk.Menu(self, tearoff=0)
        menu.add_command(label="Paste", command=(lambda :self.paste(event)))
        menu.tk_popup(event.x_root, event.y_root)
        self.bind("<ButtonRelease-3>", lambda event: None)

    #############################
    #  SELECTING
    ##############################

    def start_selection_field(self, event):
        '''
        Called when user clicks on the grid or the background of the canvas.
        Builds the dashed selection box (initially with zero area) 
        where the one corner is located at the click position, and the other
        is located at the current mouses position.
        The box will be deleted when the click is released.
        '''

        self.deselect_all()

        # Store location at which the user clicks 
        # in canvas units.
        # This will form one corner of the selection box.
        self.selection_rectangle_x_start = self.canvasx(event.x)
        self.selection_rectangle_y_start = self.canvasy(event.y)

        # Create the dashed box.
        # The other corner of the selection box is determined by the position
        # the mouse (for the moment the box has zero area).
        self.selection_rectangle = self.create_rectangle(
            self.canvasx(event.x), self.canvasy(event.y), self.canvasx(event.x), self.canvasy(event.y), 
            dash=(3, 5))

    def expand_selection_field(self, event):
        '''
        Called when user clicks+drags the mouse
        on the grid or the background of the canvas.
        Will continuously deselct all the components, then go through all
        the components and selecting those contained in the selection box.
        Wheter a component is in or out of the box is determined by the components
        box_select method which takes the coordinates of the selection rectangle as arguments.
        '''
        self.deselect_all()
        self.coords(self.selection_rectangle,
                    min(self.canvasx(event.x), self.selection_rectangle_x_start),
                    min(self.canvasy(event.y), self.selection_rectangle_y_start),
                    max(self.canvasx(event.x), self.selection_rectangle_x_start),
                    max(self.canvasy(event.y), self.selection_rectangle_y_start))
        for el in self.elements:
            el.box_select(*self.coords(self.selection_rectangle))

    def end_selection_field(self, event):
        '''
        Called when user releases a click on the grid or background of the canvas.
        Will delete the selection rectangle.
        '''
        self.delete(self.selection_rectangle)

    def deselect_all(self, event=None):
        '''
        Deselects all components on the canvas.
        '''
        for el in self.elements:
            el.deselect()

    def select_all(self, event=None):
        '''
        Selects all components on the canvas.
        '''
        for el in self.elements:
            el.force_select()

    #############################
    #  DELETING
    ##############################

    def delete_selection(self, event=None, track_changes=None):
        '''
        Deletes selected components.

        Parameters
        ----------
        track_changes:  Boolean or None
                        default is None, in which case the function track changes according
                        to the value of self.track_changes
                        True (or False), force the function to keep track (or not) of the
                        deletion of all components, allowing the user to undo (or not) this
                        operation.
        '''

        was_tracking_changes = self.track_changes

        # Delete all selected components without tracking those changes
        self.track_changes = False
        to_delete = [el for el in self.elements if el.selected]
        for el in to_delete:
            el.delete()

        # Append these changes to the history variable depending on the value 
        # of the track_changes input parameters
        if track_changes is None:  # Just follow "was_tracking_changes"
            self.track_changes = was_tracking_changes
            self.save()
        elif track_changes is False:
            self.track_changes = was_tracking_changes
        elif track_changes is True:
            self.track_changes = True
            self.save()
            self.track_changes = track_changes

        # This may have deleted components out of the visible canvas, 
        # meaning that the scrollable region should be re-configured.
        self.configure_scrollregion()

    def delete_all(self, event=None, track_changes=None):
        '''
        Selects all components, then deletes the selection

        Parameters
        ----------
        track_changes:  Boolean or None
                        default is None, in which case the function track changes according
                        to the value of self.track_changes
                        True (or False), force the function to keep track (or not) of the
                        deletion of all components, allowing the user to undo (or not) this
                        operation.
        '''
        self.select_all()
        self.delete_selection(event, track_changes)

    #############################
    #  CIRCLE utilities
    ##############################

    def create_circle(self, x, y, r):
        '''
        Creates a filled black circle located at 
        (x,y) with a radius r. 

        Parameters
        ----------
        x:  float
            horizontal position of the center of the circle in canvas units
        y:  float
            vertical position of the center of the circle in canvas units
        r:  float
            radius of the circle in canvas units

        Returns
        -------
        object_id:  integer
            object ID of the circle (actually an tkinter oval)
        '''

        x0 = x - r
        y0 = y - r
        x1 = x + r
        y1 = y + r
        return self.create_oval(x0, y0, x1, y1, fill='black')

    def update_circle(self, circle, x, y, r):
        '''
        Updates the position and size of a circle.

        Parameters
        ----------
        circle: int
            object ID of a circle (actually an tkinter oval)
        x:  float
            horizontal position of the center of the circle in canvas units
        y:  float
            vertical position of the center of the circle in canvas units
        r:  float
            radius of the circle in canvas units
        '''

        x0 = x - r
        y0 = y - r
        x1 = x + r
        y1 = y + r
        self.coords(circle, x0, y0, x1, y1)

    #############################
    #  POSITIONNING
    ##############################

    def grid_to_canvas(self, pos):
        '''
        Converts a position in grid units to canvas units

        Parameters
        ----------
        pos:    list of floats
                pos = [x_grid,y_grid] where x_grid and y_grid are given in grid units

        Returns
        -------
        [x_canvas,y_canvas]:    list of floats
                                position in canvas units
        '''

        return [self.canvas_center[0]+self.grid_unit*pos[0],
                self.canvas_center[1]+self.grid_unit*pos[1]]

    def canvas_to_grid(self, pos):
        '''
        Converts a position in canvas units to grid units

        Parameters
        ----------
        pos:    list of floats
                pos = [x_canvas,y_canvas] where x_canvas and y_canvas are given in canvas units

        Returns
        -------
        [x_grid,y_grid]:    list of floats
                            position in grid units
        '''

        return [(pos[0]-self.canvas_center[0])/self.grid_unit,
                (pos[1]-self.canvas_center[1])/self.grid_unit]
    
    def get_mouse_location(self):
        '''
        Returns the location of the mouse pointer in canvas units.

        Returns
        -------
        [x,y]:  List of floats
                corresponding to the location of the mouse pointer in canvas units.
        '''

        return [self.canvasx(self.winfo_pointerx())-self.winfo_rootx(), 
            self.canvasy(self.winfo_pointery())-self.winfo_rooty()]
    
    def center_window_on_circuit(self):
        '''
        Called when the gui is opened or a new circuit is opened, or
        upon clicking View>Re-center.

        Moves the canvas such that the top left point of the circuit 
        is situated at the top left of the canvas.
        '''

        margin = 3
        '''
        Mrgin between the circuit and the upper and left edges of the canvas.
        In grid units.
        '''

        if len(self.elements) > 0:

            # determine the coordinates of a box surrounding all
            # the circuit elements (plus a margin)
            xs = [el.x_minus for el in self.elements] + \
                [el.x_plus for el in self.elements]
            ys = [el.y_minus for el in self.elements] + \
                [el.y_plus for el in self.elements]
            box_elements = self.grid_to_canvas([min(xs)-margin, min(ys)-margin])\
                    +self.grid_to_canvas([max(xs)+margin, max(ys)+margin])

            # Set this box to be the scrollable region and
            # move to the upper left corner of that scrollable region
            self.configure(scrollregion=box_elements)
            self.xview_moveto(0)
            self.yview_moveto(0)

        # Since the canvas has moved, 
        # re-configure the scrollregion and 
        # re-draw the grid
        self.configure_scrollregion()
        self.draw_grid()

    #############################
    #  UTILITIES
    ##############################
    def message(self, text, t = 0.3, x = 5, y = 2,size = 8, weight = 'normal'):
        '''
        Displays a message on the canvas for a given amount of time.
        Called each time the circuit is saved, or to 
        inform the user he cannot zoom anymore, etc...

        Parameters
        ----------
        text:   string
                message to be displayed
        t:      float, optional
                Amount of time to display the message (in seconds).
                Default is 0.3
        x:      float, optional
                Distance of the message with respect to the left
                side of the window, in window units. Default is 5
        y:      float, optional
                Distance of the message with respect to the top
                side of the window, in window units. Default is 2
        size:   float, optional
                Fontsize to use. Default is 8
        weight: string, optional
                Font weight. Default is 'normal'
        '''

        saved_message = self.create_text(
            self.canvasx(x), 
            self.canvasy(y), 
            text=text, 
            anchor=tk.NW,
            font=Font(family='Helvetica', 
            size=size, 
            weight=weight))

        self.after(int(1000*t), lambda: self.delete(saved_message))

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

    def abort_creation(self,event = None, rerun_command = True):
        self.canvas.in_creation = None
        self.canvas.set_keyboard_shortcuts_element_creation()
        if event.type == tk.EventType.KeyPress and rerun_command:
            self.canvas.event_generate(event.char)
        del self


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

    def release_motion_paste_single(self, event):
        self.release_motion_paste(event)
        self.canvas.bind('<Left>', lambda event: None)
        self.canvas.bind('<Right>', lambda event: None)
        self.canvas.bind('<Up>', lambda event: None)
        self.canvas.bind('<Down>', lambda event: None)

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
        
    def add_nodes(self, to = 'all wires', minus = True, plus = True):

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
            if w.x_minus == w.x_plus == self.x_minus and minus:
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

            elif w.y_minus == w.y_plus == self.y_minus and minus:
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
            elif w.x_minus == w.x_plus == self.x_plus and plus:
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

            elif w.y_minus == w.y_plus == self.y_plus and plus:
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

    def add_or_replace_node_dots(self, plus = True, minus = True):
        gu = self.canvas.grid_unit
        canvas_coords_minus = self.grid_to_canvas([self.x_minus,self.y_minus])

        if minus:
            if self.dot_minus is None:
                self.dot_minus = self.canvas.create_circle(
                    *canvas_coords_minus, gu*node_dot_radius)
            else:
                self.canvas.update_circle(self.dot_minus,
                    *canvas_coords_minus, gu*node_dot_radius)
            self.canvas.tag_raise(self.dot_minus)

        if plus:
            canvas_coords_plus = self.grid_to_canvas([self.x_plus,self.y_plus])
            if self.dot_plus is None:
                self.dot_plus = self.canvas.create_circle(
                    *canvas_coords_plus, gu*node_dot_radius)
            else:
                self.canvas.update_circle(self.dot_plus,
                    *canvas_coords_plus, gu*node_dot_radius)

            self.canvas.tag_raise(self.dot_plus)

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
        menu.tk_popup(event.x_root, event.y_root)
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
        self.canvas.set_keyboard_shortcuts_element_creation()
        self.canvas.in_creation = None
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
        self.canvas.bind("<Escape>", lambda event: self.abort_creation(event, rerun_command = False))
        self.canvas.bind('r', self.abort_creation)
        self.canvas.bind('l', self.abort_creation)
        self.canvas.bind('c', self.abort_creation)
        self.canvas.bind('j', self.abort_creation)
        self.canvas.bind('w', self.abort_creation)
        self.canvas.bind('g', self.abort_creation)

    def auto_place(self, auto_place_info):
        self.create()
        self.add_nodes()

    def start_line(self, event):
        self.x_minus, self.y_minus = self.init_minus_snap_to_grid(event)

        gu = self.canvas.grid_unit
        self.dot_minus = self.canvas.create_circle(
            *(self.canvas.grid_to_canvas([self.x_minus,self.y_minus])), gu*node_dot_radius)

        self.canvas.bind("<Motion>", self.show_line)
        self.canvas.bind("<Button-1>", self.end_line)
        self.canvas.in_creation = self

    def abort_creation(self, event=None, rerun_command = True):
        self.canvas.bind("<Button-1>", lambda event: None)
        self.canvas.bind("<Escape>", lambda event: None)
        self.canvas.config(cursor='arrow')
        if self.dot_minus is not None:
            # First node has been created
            self.canvas.bind("<Motion>", lambda event: None)
            self.canvas.delete('temp')
            self.canvas.delete(self.dot_minus)
        super(W,self).abort_creation(event,rerun_command)

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

    def init_adapt_to_grid_unit(self, event):
        self.show_line(event)
        self.add_or_replace_node_dots(plus = False)

    def show_line(self, event):
        self.canvas.delete("temp")

        xm, ym = self.grid_to_canvas([self.x_minus, self.y_minus])
        xp = self.canvas.canvasx(event.x)
        yp = self.canvas.canvasy(event.y)

        if abs(xm-xp) > abs(ym-yp):
            # Horizontal line
            self.line = self.canvas.create_line(xm, ym, xp, ym, tags='temp',
            width=lw*self.canvas.grid_unit,
            fill = light_black)
        else:
            # Vertical line
            self.line = self.canvas.create_line(xm, ym, xm, yp, tags='temp',
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
        self.canvas.in_creation = self

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

    def import_image(self):
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
        size = int(self.canvas.grid_unit*(1-1*node_dot_radius))
        img = img.resize((size, int(size/2)))
        img = img.rotate(angle,expand = True)
        self.tk_image = ImageTk.PhotoImage(img)

    def create(self):
        self.add_or_replace_node_dots()
        x, y, angle = self.pos
        self.import_image()
        self.image = self.canvas.create_image(
            *self.grid_to_canvas([x, y]), image=self.tk_image)
        self.add_or_replace_label()
        self.canvas.elements.append(self)
        self.set_allstate_bindings()

    def update_graphic(self):
        self.import_image()
        self.canvas.itemconfig(self.image, image=self.tk_image)

    def adapt_to_grid_unit(self):
        self.update_graphic()
        self.canvas.coords(self.image, *self.grid_to_canvas(self.pos[:2]))
        self.add_or_replace_label()
        self.add_or_replace_node_dots()
        
    def init_adapt_to_grid_unit(self,event):
        self.update_graphic()

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

        self.init_angle = angle
        self.import_image()
        self.image = self.canvas.create_image(
            self.canvas.canvasx(event.x), self.canvas.canvasy(event.y), image=self.tk_image)

        self.canvas.bind("<Button-1>", self.init_release)
        self.canvas.bind('<Motion>', self.init_on_motion)
        self.canvas.bind('<Escape>', lambda event: self.abort_creation(event, rerun_command = False))
        self.canvas.bind('r', self.abort_creation)
        self.canvas.bind('l', self.abort_creation)
        self.canvas.bind('c', self.abort_creation)
        self.canvas.bind('j', self.abort_creation)
        self.canvas.bind('w', self.abort_creation)
        self.canvas.bind('g', self.abort_creation)
        self.canvas.bind(
            '<Left>', lambda event: self.init_create_component(event, angle=WEST))
        self.canvas.bind(
            '<Right>', lambda event: self.init_create_component(event, angle=EAST))
        self.canvas.bind(
            '<Up>', lambda event: self.init_create_component(event, angle=NORTH))
        self.canvas.bind(
            '<Down>', lambda event: self.init_create_component(event, angle=SOUTH))

    def unset_initialization_bindings(self):
        self.canvas.bind("<Button-1>", lambda event: None)
        self.canvas.bind('<Motion>', lambda event: None)
        self.canvas.bind('<Left>', lambda event: None)
        self.canvas.bind('<Right>', lambda event: None)
        self.canvas.bind('<Up>', lambda event: None)
        self.canvas.bind('<Down>', lambda event: None)
        self.canvas.bind('<Escape>', lambda event: None)

    def abort_creation(self, event=None, rerun_command = True):
        self.unset_initialization_bindings()
        self.canvas.delete(self.image)
        self.canvas.delete(self.dot_minus)
        self.canvas.delete(self.dot_plus)
        super(Component,self).abort_creation(event, rerun_command)

    def init_release(self, event):
        self.unset_initialization_bindings()
        self.canvas.track_changes = False
        self.snap_to_grid()
        self.canvas.set_keyboard_shortcuts_element_creation()
        self.request_value_label()
        self.canvas.in_creation = None
        if self.prop[0] is None and self.prop[1] is None:
            self.abort_creation(rerun_command = False)
            self.canvas.track_changes = True
            return
        self.add_or_replace_label()
        self.canvas.elements.append(self)
        self.set_allstate_bindings()
        self.canvas.track_changes = True
        self.canvas.save()

        # If mouse is still on top of component,
        # act as if one had hovered on top of component

        x,y = self.canvas.get_mouse_location()

        # bounding box of image
        xm,ym,xp,yp = self.canvas.bbox(self.image)

        if xm<x<xp and ym<y<yp:
            self.hover_enter(event)

    def open_right_click_menu(self, event):
        menu = tk.Menu(self.canvas, tearoff=0)
        menu.add_command(label="Edit", command=self.modify_values)
        menu.add_command(label="Rotate", command=self.rotate)
        menu.add_command(label="Delete", command=self.canvas.delete_selection)
        menu.add_separator()
        menu.add_command(label="Copy", command=self.canvas.copy_selection)
        menu.add_command(label="Cut", command=self.canvas.cut_selection)
        menu.tk_popup(event.x_root, event.y_root)
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
        menu.tk_popup(event.x_root, event.y_root)
        self.canvas.bind("<ButtonRelease-3>", lambda event: None)

    def add_or_replace_label(self):
        pass

    def add_or_replace_node_dots(self):
        super(G,self).add_or_replace_node_dots(minus = False)

    def add_nodes(self, to = 'all wires', minus = False, plus = True):   
        super(G,self).add_nodes(to = to, minus = minus, plus = plus)       

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

class GuiWindow(ttk.Frame):
    """
    GuiWindow inherits from the tkinter Frame class.
    A Frame is a rectangular region on the screen.
    This Frame just plays the role of placeholder for another 
    Frame defined in the CircuitEditor which will
    host the canvas, menubar, scrollbars.

    In this class we manage the launching of the GUI
    and properties pertaining to the opened window
    such as the titlebar and initial window size
    
    Parameters
    ----------
    netlist_filename:   string
                        path to the file used to save the network constructed
                        in the GUI
    """

    def __init__(self, netlist_filename):

        # Initialize the frame, inside the root window (tk.Tk())
        ttk.Frame.__init__(self, master=tk.Tk())

        # Set the name to appear in the title bar
        self.master.title('Circuit Editor')

        # Set the initial size of the window in pixels
        self.master.geometry('800x600')

        # Load the logo to the title bar
        try:
            self.master.iconbitmap(r'C:\ProgramData\Anaconda3\Lib\site-packages\Qcircuits\artwork\logo.ico')
        except Exception as e:
            # Anticipating possible non-Windows related issues
            print("There has been an error loading the applications icon:\n"+str(e))

        # Make the fram a 1x1 expandable grid
        self.master.rowconfigure(0, weight=1)
        self.master.columnconfigure(0, weight=1)

        # Populate that grid with the circuit editor
        self.canvas = CircuitEditor(
            self.master, netlist_filename=netlist_filename, grid_unit=60)
        self.mainloop()

if __name__ == '__main__':
    GuiWindow('./src/test.txt')
