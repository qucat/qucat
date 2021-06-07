import sys
import os
sys.path.append(os.path.join(os.path.join(os.path.dirname(os.path.dirname(__file__)),'src')))
from utils import run_notebook, TestCaseAppended
import unittest
from shutil import rmtree

class TestTutorials(TestCaseAppended):
    def run_tutorial(self,notebook_name):
        var_dict = run_notebook(notebook_name)
        self.assertFalse(isinstance(var_dict,str),msg=var_dict)
        return var_dict


class TestAll(TestTutorials):

    def test_basics(self):
        var_dict = self.run_tutorial('basics.ipynb')
    def test_filter(self):
        var_dict = self.run_tutorial('filter_design.ipynb')
    def test_mmusc(self):
        var_dict = self.run_tutorial('MMUSC.ipynb')
    def test_OM(self):
        var_dict = self.run_tutorial('optomechanics.ipynb')
    def test_TC(self):
        var_dict = self.run_tutorial('tuneable_coupler.ipynb')
    def test_SNAIL(self):
        var_dict = self.run_tutorial('snail.ipynb')
        
if __name__ == "__main__":
    unittest.main()