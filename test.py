import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__),'src'))

import numpy as np

from core import Network,GUI,J,L,C,R
from time import time
import matplotlib.pyplot as plt
c = GUI('test.txt', edit = True,plot = False)
c.f_k_A_chi(pretty_print=True)
