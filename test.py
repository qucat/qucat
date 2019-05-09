import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__),'src'))

from core import Network,GUI,J,L,C,R
# Cj = 100e-15
# Lj = 10e-9
# junction = J(0,1,'Lj')
# circuit = Network([
#     C(0,1,Cj),
#     junction
#     ])
# circuit.f_k_A_chi(Lj=1)
# junction.zpf(mode=0,quantity = 'flux')
# H = circuit.hamiltonian(modes = [0],taylor = 4,excitations = [50])
# print(H)
circuit = GUI(filename = 'test.txt',edit=False,plot=False)
circuit.show_normal_mode(0,quantity='current',L_J=1e-8)
# circuit.f_k_A_chi()
# print(circuit.resistors[0].phasor(0,'voltage'))
# circuit.hamiltonian(L_J = 1e-9,modes=[0],excitations=[5],return_ops=True,taylor=4)
# circuit.eigenfrequencies(L_J = np.linspace(1e-9,2e-9,4))
# circuit.f_k_A_chi(L_J = np.linspace(1e-9,2e-9,4))
# print(circuit.Y)
# print(sp.together(circuit.Y))
# print(circuit.eigenfrequencies())


