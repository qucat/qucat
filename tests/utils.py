import traceback
import nbformat  # to convert the notebook into an array of cells
from os.path import join, dirname
from os import devnull
import sys

sys.path.append(join(join(dirname(dirname(__file__)), "src")))
import unittest
import core
import numpy as np

# Run plt.ion() to avoid hanging on plt.show() calls
import matplotlib.pyplot as plt

plt.ion()


def cutoff_digits(f, digits):
    float_format = "%%.%de" % digits
    return float(float_format % (np.real(f))) + 1j * float(float_format % (np.imag(f)))


class TestCaseAppended(unittest.TestCase):
    def assertRelativelyClose(self, a, b, digits=6):
        a = cutoff_digits(a, digits)
        b = cutoff_digits(b, digits)
        self.assertEqual(a, b)

    def assertArrayRelativelyClose(self, a, b, digits=6):
        a = np.array(a)
        b = np.array(b)
        self.assertTrue(
            a.shape == b.shape,
            msg=f"Arrays do not have the same dimension {a.shape}!={b.shape}",
        )
        for index, _ in np.ndenumerate(a):
            a_comp = cutoff_digits(a[index], digits)
            b_comp = cutoff_digits(b[index], digits)
            self.assertEqual(
                a_comp,
                b_comp,
                msg=f"Components with index {index} do not match {a_comp}!={b_comp}",
            )

    def open_gui_file(self, filename, edit=False, print_network=False, plot=False):
        return core.GUI(
            join(dirname(__file__), "gui_testing_files", filename),
            edit=edit,
            print_network=print_network,
            plot=plot,
        )


def parse_cell(code):
    lines = code.split("\n")

    valid_lines = ""
    for line in lines:
        if len(line.replace(" ", "")) > 0:
            if not line.replace(" ", "")[0] == "%" and not "plt.show(" in line:
                if "GUI(" in line:
                    valid_lines += "def GUI(*args,**kwargs):\n"
                    valid_lines += "    from qucat import GUI\n"
                    valid_lines += "    return GUI(args[0],edit=False)\n"
                valid_lines += line + "\n"

        valid_lines = valid_lines.replace("qucat.", "").replace(
            "from qucat ", "from core "
        )
    return valid_lines


def run_notebook(notebook_name):
    import sys
    import os

    sys.path.append(join(join(dirname(dirname(__file__)), "src")))
    tutorials_folder = join(
        dirname(dirname(__file__)), "docs_src", "source", "tutorials"
    )
    os.chdir(tutorials_folder)
    nb = nbformat.read(join(tutorials_folder, notebook_name), 4)
    sys.stdout = open(devnull, "w")
    sys.stderr = open(devnull, "w")
    for c in nb["cells"]:
        if c["cell_type"] == "code":
            to_run = parse_cell(c["source"])
            try:
                exec(to_run, globals())
            except Exception as e:
                tb = "\nThere was an error executing the following notebook code:\n"
                tb += "==============================\n"
                n = 1
                for line in to_run.splitlines():
                    tb += "% 4d %s \n" % (n, line)
                    n += 1
                tb += "==============================\n"
                tb += traceback.format_exc()
                sys.stdout = sys.__stdout__
                return tb
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
    return dict(globals())
