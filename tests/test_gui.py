import unittest
import os
from Qcircuits.src._gui import GuiWindow
import inspect

def get_netlist_filename():
    curframe = inspect.currentframe()
    calframe = inspect.getouterframes(curframe, 2)
    calling_function_name = calframe[2][3]
    filename = os.path.join(\
        os.path.dirname(__file__),\
        ".temp",\
        calling_function_name+\
        "_netlist.txt")
    return filename

def write_netlist_file(contents):
    filename = get_netlist_filename()
    with open(filename,'w') as netlist_file:
        netlist_file.write(contents)
    return filename

def read_netlist_file():
    with open(get_netlist_filename(),'r') as netlist_file:
        contents = netlist_file.read()
    return contents

class TestOpening(unittest.TestCase):

    def test_if_opening_blank_test_throws_error(self):
        filename = write_netlist_file('')
        gui = GuiWindow(filename, unittesting = True)
        gui.update()
        gui.destroy()
        self.assertEqual('',read_netlist_file())