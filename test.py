import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__),'src'))

import numpy as np

from core import Network,GUI,J,L,C,R
from time import time
import matplotlib.pyplot as plt
c = GUI('test.txt', edit = False,plot = False)
# kwargs = {'Cc2':1e-25,'C1':165e-15,'L2':10e-9,'C2':100e-15,'L1':2.7e-9,'R':50,'Cc1':1e-16}
# c.f_k_A_chi(pretty_print=True,**kwargs)
C_list = np.logspace(-15,-25,101)
to_plot = []
for i,C in enumerate(C_list):
    kwargs = {'Cc2':C,'C1':165e-15,'L2':10e-9,'C2':100e-15,'L1':2.7e-9,'R':50,'Cc1':1e-16}
    freqs,_,_,_ = c.f_k_A_chi(pretty_print=False,**kwargs)
    # to_plot.append(c.anharmonicities(**kwargs)[1])
    to_plot.append([c.components['L1']._flux_zpf_r(c.zeta[1],**kwargs),
        c.components['L2']._flux_zpf_r(c.zeta[1],**kwargs)])
print (to_plot)
plt.loglog(C_list,to_plot)
plt.show()