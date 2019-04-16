import unittest
import os
import shutil
from Qcircuits.src._gui import GuiWindow
import inspect

class ManualTesting(unittest.TestCase):

    def get_netlist_filename(self):
        curframe = inspect.currentframe()
        calframe = inspect.getouterframes(curframe, 2)
        calling_function_name = calframe[2][3]
        filename = os.path.join(\
            os.path.dirname(__file__),\
            ".gui_testing_files",\
            calling_function_name+\
            "_netlist.txt")
        return filename

    def write_netlist_file(self,contents):
        filename = self.get_netlist_filename()
        with open(filename,'w') as netlist_file:
            netlist_file.write(contents)
        return filename

    def read_netlist_file(self):
        with open(self.get_netlist_filename(),'r') as netlist_file:
            contents = netlist_file.read()
        return contents


class GuiTesting(unittest.TestCase):

    def launch_gui_testing(self,exclude = None):
        folder = self.get_folder_name()

        if not self.already_built(folder):
            self.gui_build_test(folder,exclude)

        init = os.path.join(folder,'initial_netlist.txt')
        events = os.path.join(folder,'events.txt')
        self.run_events(folder,init,events)

        with open(os.path.join(folder,'final_netlist.txt'),'r') as f:
            final_expected = f.read()
        with open(os.path.join(folder,'final_after_events_netlist.txt'),'r') as f:
            final_after_events = f.read()
            
        self.assertEqual(final_expected,final_after_events)

    def get_folder_name(self):
        curframe = inspect.currentframe()
        calframe = inspect.getouterframes(curframe, 2)
        calling_function_name = calframe[2][3]
                
        folder = os.path.join(\
            os.path.dirname(__file__),\
            ".gui_testing_files",\
            calling_function_name)
        return folder

    def run_events(self,folder,  init, events):
        init = os.path.join(folder,'initial_netlist.txt')
        final_after_events = os.path.join(folder,'final_after_events_netlist.txt')
        events = os.path.join(folder,'events.txt')

        shutil.copyfile(init,final_after_events)
        gui = GuiWindow(final_after_events, _unittesting = True)
        with open(events,'r') as f:
            lines = f.readlines()
            for l in lines:
                if 'Motion' in l:
                    l += ',warp=True'
                exec('gui.canvas.event_generate('+l+')', globals(), locals())
        gui.destroy()



    def gui_build_test(self,folder,exclude):
        init = os.path.join(folder,'initial_netlist.txt')
        final = os.path.join(folder,'final_netlist.txt')
        events = os.path.join(folder,'events.txt')

        print("Build initial circuit")
        GuiWindow(init)
        shutil.copyfile(init,final)
        print("Build final circuit")
        GuiWindow(final, _track_events_to = events)

        self.remove_exluded_events(folder,exclude)

    def remove_exluded_events(self,folder,exclude):
        if exclude is not None:
            events = os.path.join(folder,'events.txt')
            temp = os.path.join(folder,'temp.txt')
            shutil.copyfile(events,temp)
            with open(events,'w') as to_write:
                with open(temp,'r') as to_read:
                    lines = to_read.readlines()
                    for l in lines:
                        wr = True
                        for e in exclude:
                            if e in l:
                                wr = False
                        if wr:
                            to_write.write(l)
            os.remove(temp)


    def already_built(self,folder):
        try:
            os.mkdir(folder)
            return False
        except FileExistsError:
            pass
        
        for filename in [
            'initial_netlist.txt',
            'final_netlist.txt',
            'events.txt']:
            try:
                with open(os.path.join(folder,filename)):
                    pass
            except FileNotFoundError:
                return False 
        return True


class TestOpening(ManualTesting):

    def test_if_opening_blank_test_throws_error(self):
        filename = self.write_netlist_file('')
        gui = GuiWindow(filename, _unittesting = True)
        gui.destroy()
        self.assertEqual('',self.read_netlist_file())

class TestMovingComponentsAround(GuiTesting):

    def test_moving_capacitor_horizontally(self):
        self.launch_gui_testing()

    def test_rotating_capacitor(self):
        self.launch_gui_testing()

if __name__ == "__main__":
    unittest.main()
    # unittest.main(defaultTest='TestMovingComponentsAround')
    # unittest.main(defaultTest='TestMovingComponentsAround.test_rotating_capacitor')