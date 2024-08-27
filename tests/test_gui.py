import unittest
import os
import sys
import shutil

# We have to do this to avoid calling the package __init__.py
# which calls matplotlib.pyplot, causing issues on MAC OSX
# see https://stackoverflow.com/questions/32019556/matplotlib-crashing-tkinter-application
sys.path.append(
    os.path.join(
        os.path.join(os.path.dirname(os.path.dirname(__file__)), "src")
    )
)
from _gui import GuiWindow
import inspect
import tkinter as tk


class GuiTestingHandler(unittest.TestCase):
    pass


class ManualTesting(GuiTestingHandler):
    """
    Example code:

    def test_if_opening_blank_test_throws_error(self):
        filename = self.write_netlist_file('')
        self.gui = GuiWindow(filename, _unittesting = True)
        self.gui.master.destroy()
        self.assertEqual('',self.read_netlist_file())
    """

    def get_netlist_filename(self):
        curframe = inspect.currentframe()
        calframe = inspect.getouterframes(curframe, 2)
        calling_function_name = calframe[2][3]
        filename = os.path.join(
            os.path.dirname(__file__),
            "gui_testing_files",
            calling_function_name + "_netlist.txt",
        )
        return filename

    def write_netlist_file(self, contents):
        filename = self.get_netlist_filename()
        with open(filename, "w") as netlist_file:
            netlist_file.write(contents)
        return filename

    def read_netlist_file(self):
        with open(self.get_netlist_filename(), "r") as netlist_file:
            contents = netlist_file.read()
        return contents


class AutomaticTesting(GuiTestingHandler):
    def launch_gui_testing(
        self,
        exclude=None,
        force_build=False,
        run_slower=False,
        os_type="windows",
    ):
        self.os_type = os_type
        self.force_build = force_build
        self.set_folder_name()
        self.set_file_names()
        self.exclude = exclude
        self.run_slower = run_slower

        if not self.already_built():
            self.gui_build_test()

        self.run_events()

        with open(self.final_expected, "r") as f:
            final_expected = f.read()
        with open(self.final_after_events, "r") as f:
            final_after_events = f.read()

        self.assertEqual(final_expected, final_after_events)

    def set_file_names(self):
        self.init = os.path.join(self.folder, "initial_netlist.txt")
        self.events = os.path.join(self.folder, "events.txt")
        self.final_expected = os.path.join(self.folder, "final_netlist.txt")
        self.final_after_events = os.path.join(
            self.folder, "final_after_events_netlist.txt"
        )

        if self.force_build:
            for filename in [
                self.init,
                self.final_expected,
                self.events,
                self.final_after_events,
            ]:
                try:
                    os.remove(filename)
                except FileNotFoundError:
                    pass

    def set_folder_name(self):
        curframe = inspect.currentframe()
        calframe = inspect.getouterframes(curframe, 2)
        calling_function_name = calframe[2][3]

        self.folder = os.path.join(
            os.path.dirname(__file__),
            "gui_testing_files",
            calling_function_name,
        )

    def run_events(self):

        shutil.copyfile(self.init, self.final_after_events)
        self.gui = GuiWindow(
            self.final_after_events, _unittesting=True, _os_type=self.os_type
        )

        with open(self.events, "r") as f:
            lines = f.readlines()
            for l in lines:
                self.gui.canvas.focus_force()
                if l[0] == "#":
                    # Scrollbar event
                    exec(
                        "self.gui.canvas.scroll_" + l[1] + "(" + l[3:] + ")",
                        globals(),
                        locals(),
                    )
                # elif '<Motion>' in l:
                #      exec('self.gui.canvas.event_generate('+l+',warp=True)', globals(), locals())
                else:
                    exec(
                        "self.gui.canvas.event_generate(" + l + ")",
                        globals(),
                        locals(),
                    )

                if self.run_slower:
                    self.gui.update()
                    self.gui.canvas.after(10)

        if self.run_slower:
            self.gui.master.mainloop()
        else:
            self.gui.master.destroy()

    def gui_build_test(self):

        print("Build initial circuit")
        GuiWindow(self.init)
        shutil.copyfile(self.init, self.final_expected)
        print("Build final circuit")
        GuiWindow(self.final_expected, _track_events_to=self.events)

        self.remove_exluded_events()

    def remove_exluded_events(self):
        if self.exclude is not None:
            temp = os.path.join(self.folder, "temp.txt")
            shutil.copyfile(self.events, temp)
            with open(self.events, "w") as to_write:
                with open(temp, "r") as to_read:
                    lines = to_read.readlines()
                    for l in lines:
                        wr = True
                        for e in self.exclude:
                            if e in l:
                                wr = False
                        if wr:
                            to_write.write(l)
            os.remove(temp)

    def already_built(self):
        try:
            os.mkdir(self.folder)
            return False
        except FileExistsError:
            pass

        for filename in [self.init, self.final_expected, self.events]:
            try:
                with open(filename, "r"):
                    pass
            except FileNotFoundError:
                return False

        return True


class TestComponentCreation(AutomaticTesting):
    # def test_building_resistor(self):
    #     self.launch_gui_testing()
    def test_building_wire(self):
        self.launch_gui_testing()

    def test_building_overlapping_and_intersecting_wires(self):
        self.launch_gui_testing()

    def test_building_capacitor(self):
        self.launch_gui_testing()

    def test_building_junction(self):
        self.launch_gui_testing()

    def test_building_ground(self):
        self.launch_gui_testing()

    def test_building_rotated_ground(self):
        self.launch_gui_testing()

    def test_building_inductor(self):
        self.launch_gui_testing()

    def test_building_transmon(self):
        self.launch_gui_testing()

    def test_cancel_wire_build_before_first_node(self):
        self.launch_gui_testing()

    def test_cancel_wire_build_after_first_node(self):
        self.launch_gui_testing()

    def test_cancel_junction_build(self):
        self.launch_gui_testing()


class TestDeleting(AutomaticTesting):
    def test_deleting_an_RLC_on_mac(self):
        self.launch_gui_testing(os_type="mac")


class TestOpening(AutomaticTesting):
    def test_if_opening_blank_test_throws_error(self):
        self.launch_gui_testing()

    def test_opening_complicated_wire_network(self):
        self.launch_gui_testing()

    def test_opening_complicated_circuit(self):
        self.launch_gui_testing()

    def test_opening_longfloat(self):
        self.launch_gui_testing()


class TestCutCopyPaste(AutomaticTesting):
    def test_copy_paste__copying_nothing(self):
        self.launch_gui_testing(run_slower=False)

    def test_cut_paste__cut_paste_inductor(self):
        self.launch_gui_testing()

    def test_copy_paste__copy_paste_capacitor(self):
        self.launch_gui_testing()

    def test_copy_paste__paste_capacitor_to_intersect_wire(self):
        self.launch_gui_testing()

    def test_cut_paste__box_select_cut_paste_random_complicated_circuit(self):
        self.launch_gui_testing(run_slower=False)

    def test_cut_paste__select_all_cut_paste_random_complicated_circuit(self):
        self.launch_gui_testing()

    def test_copy_paste__select_all_copy_paste_random_complicated_circuit(
        self,
    ):
        self.launch_gui_testing()


class TestUndoRedo(AutomaticTesting):
    def test_undo_deletion(self):
        self.launch_gui_testing(run_slower=False)

    def test_undo_movement_of_capacitor(self):
        self.launch_gui_testing(run_slower=False)

    def test_undo_deletion_redo(self):
        self.launch_gui_testing(run_slower=False)


class TestMovingComponentsAround(AutomaticTesting):
    def test_moving_capacitor_horizontally(self):
        self.launch_gui_testing(run_slower=False)

    def test_rotating_capacitor(self):
        self.launch_gui_testing(force_build=False, run_slower=False)

    def test_rotating_ground(self):
        self.launch_gui_testing(force_build=False, run_slower=False)

    def test_rotating_RLCJG_right_click(self):
        self.launch_gui_testing(force_build=False, run_slower=False)

    def test_rotating_RLCJG_draggin_arrows(self):
        self.launch_gui_testing(force_build=False, run_slower=False)

    def test_rotating_nothing_using_Alt_R(self):
        self.launch_gui_testing(force_build=False, run_slower=False)

    def test_moving_capacitor_twice(self):
        self.launch_gui_testing()

    def test_moving_parallel_RLCJG(self):
        self.launch_gui_testing(run_slower=False)

    def test_move_ground_with_rotation(self):
        self.launch_gui_testing(run_slower=False)

    def test_moving_two_elements_and_trying_to_rotate(self):
        self.launch_gui_testing(run_slower=False)

    def test_moving_an_unselected_object_whilst_another_one_is_selected(self):
        self.launch_gui_testing(run_slower=False)


class TestZooming(AutomaticTesting):
    def test_zooming_in_out_then_moving_capacitor(self):
        self.launch_gui_testing()


class TestScrolling(AutomaticTesting):
    def test_scrolling_with_scrollbars_left_and_right_then_moving_capacitor(
        self,
    ):
        self.launch_gui_testing()

    def test_scrolling_with_mouse_left_and_right_then_moving_junction(self):
        self.launch_gui_testing(force_build=False, run_slower=False)


class TestSelection(AutomaticTesting):
    def test_selection__multiple_shift_box_selects(self):
        self.launch_gui_testing()

    def test_selection__shift_box_select_over_components_and_back(self):
        self.launch_gui_testing(run_slower=False)

    def test_selection__ctrl_shift_click_select(self):
        self.launch_gui_testing(run_slower=False, force_build=False)


if __name__ == "__main__":
    unittest.main()
