import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__),'src'))

import numpy as np

from core import Network,GUI,J,L,C,R
from time import time
# circuit = GUI(filename = 'test.txt',edit=True,plot=False)
# circuit.show_normal_mode(0,quantity='flux')
f0 = 4.603e9
Z0 = 50
Ej = 36.3e9/2
Cc = 40.3e-15
Cj = 5.13e-15

w0 = f0*2.*np.pi
C0 = np.pi/4/w0/Z0
L0 = 4*Z0/np.pi/w0

mmusc_circuits = []
dissipationless_time = []
dissipative_time = []
bare_mode_frequency = []
for N in range(9):
    netlist = [
        J(0,1,Ej,use_E=True),
        C(0,1,Cj),
        C(1,2%(2+N),Cc)]
    for m in range(N):
        node_minus = 2+m
        node_plus = (2+m+1)%(2+N)
        Lm = L0/(2*m+1)**2
        netlist += [
            L(node_minus,node_plus,Lm),
            C(node_minus,node_plus,C0)]
    try:
        ts = time()
        circuit = Network(netlist)
        circuit.f_k_A_chi()
        te = time()
        print(N)
        print(time() - ts)
    
        ts = time()
        circuit = Network(netlist+[R(2+m,(2+m+1)%(2+N),1e6) for m in range(N)])
        circuit.f_k_A_chi()
        print("with dissipation: %f"%(time() - ts))
    except Exception:
        pass

