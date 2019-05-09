import sys
import os
sys.path.append(os.path.join(os.path.join(os.path.dirname(os.path.dirname(__file__)),'src')))
from utils_test_tutorials import run_notebook
from test_core import TestCaseAppended
import unittest
from shutil import rmtree

class TestTutorials(TestCaseAppended):
    def run_tutorial(self,notebook_name):
        var_dict = run_notebook(notebook_name)
        self.assertFalse(isinstance(var_dict,str),msg=var_dict)
        try:
            rmtree('circuits')
        except FileNotFoundError:
            # no folder was created
            pass
        return var_dict


class TestBasics(TestTutorials):

    def test_transmon_LC_GUI(self):
        var_dict = self.run_tutorial('transmon_LC_GUI.ipynb')

    def test_transmon_LC_programmatically(self):
        var_dict = self.run_tutorial('transmon_LC_programmatically.ipynb')
        
if __name__ == "__main__":
    unittest.main()