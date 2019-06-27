import sys
import os
sys.path.append(os.path.join(os.path.join(os.path.dirname(os.path.dirname(__file__)),'src')))
import unittest
import _utility
from math import isclose
import numpy as np
from utils import TestCaseAppended

class All(TestCaseAppended):
    def test_dfridr(self):
        self.assertEqual(_utility.dfridr(lambda x:x**2, 0,0.1)[0],0)
        self.assertEqual(_utility.dfridr(lambda x:x**2+x, 0,0.1)[0],1)
        self.assertEqual(_utility.dfridr(lambda x:x**2+x, 0,10)[0],1)
        self.assertRelativelyClose(_utility.dfridr(lambda x:np.exp(x), 30,0.1)[0],np.exp(30))
        self.assertRelativelyClose(_utility.dfridr(lambda x:np.exp(x), 30,1)[0],np.exp(30))

        
    def test_ridders_derivative(self):
        self.assertRelativelyClose(_utility.ridders_derivative(lambda x:1j*(x**2+x), 1),3j)
        self.assertRelativelyClose(_utility.ridders_derivative(lambda x:x**2+x, 1),3)
        self.assertRelativelyClose(_utility.ridders_derivative(lambda x:np.exp(x), 30),np.exp(30))

if __name__ == "__main__":
    unittest.main()