import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__),'src'))

import numpy as np

from core import Network,GUI,J,L,C,R
from time import time
import matplotlib.pyplot as plt
c = GUI('test.txt', edit = False,plot = False)
# C_list = np.logspace(-15,-30,16)
# plt.loglog(C_list,circuit.anharmonicities(C=C_list)[1])
# plt.show()
# circuit = GUI(filename = 'test.txt',edit=False,plot=False)
kwargs = {'Cc2':1.11e-18,'C1':165e-15,'L2':10e-9,'C2':100e-15,'L1':2.7e-9,'R':50,'Cc1':1e-16}
c.f_k_A_chi(pretty_print=True,**kwargs)

# C_list = np.logspace(-15,-20,11)
# c.root_relative_tolerance = 1e-13
# to_plot = []
# for C in C_list:
#     kwargs = {'Cc2':C,'C1':165e-15,'L2':10e-9,'C2':100e-15,'L1':2.7e-9,'R':50,'Cc1':1e-16}
#     to_plot.append(c.anharmonicities(**kwargs)[1])
# plt.loglog(C_list,to_plot)
# plt.show()

# circuit.show_normal_mode(0,quantity='flux')
# f0 = 4.603e9
# Z0 = 50
# Ej = 36.3e9/2
# Cc = 40.3e-15
# Cj = 5.13e-15

# w0 = f0*2.*np.pi
# C0 = np.pi/4/w0/Z0
# L0 = 4*Z0/np.pi/w0

# mmusc_circuits = []
# dissipationless_time = []
# dissipative_time = []
# bare_mode_frequency = []
# for N in range(8):
#     netlist = [
#         J(0,1,Ej,use_E=True),
#         C(0,1,Cj),
#         C(1,2%(2+N),Cc)]
#     for m in range(N):
#         node_minus = 2+m
#         node_plus = (2+m+1)%(2+N)
#         Lm = L0/(2*m+1)**2
#         netlist += [
#             L(node_minus,node_plus,Lm),
#             C(node_minus,node_plus,C0)]
#     try:
#         ts = time()
#         circuit = Network(netlist)
#         circuit.f_k_A_chi()
#         te = time()
#         print(N)
#         print(time() - ts)
    
#         ts = time()
#         circuit = Network(netlist+[R(2+m,(2+m+1)%(2+N),1e6) for m in range(N)])
#         circuit.f_k_A_chi()
#         print("with dissipation: %f"%(time() - ts))
#     except Exception:
#         pass

