import sympy as sp
import numpy as np
from sympy.core.mul import Mul,Pow,Add

class Circuit(object):
    """docstring for Circuit"""
    def __init__(self):
        pass
    def admittance(self):
        pass
    def remove_resistances(self):
        pass

    def __add__(self, other_circuit):
        return Series(self, other_circuit)

    def __or__(self, other_circuit):
        return Parallel(self, other_circuit)



class Connection(Circuit):
    """docstring for Connection"""
    def __init__(self, component_left,component_right):
        super(Connection, self).__init__()
        self.component_left = component_left
        self.component_right = component_right
    def remove_resistances(self):
        if type(self.component_left)==R:
            return self.component_right
        elif type(circuit.component_right)==R:
            return self.component_left
        else:
            return Connection(self.component_left,
                        self.component_right)

class Series(Connection):
    """docstring for Series"""
    def __init__(self, component_left,component_right):
        super(Series, self).__init__(component_left,component_right)

    def admittance(self):
        return 1/Add(1/self.component_left.admittance(),
            1/self.component_right.admittance())

class Parallel(Connection):
    """docstring for Parallel"""
    def __init__(self, component_left,component_right):
        super(Parallel, self).__init__(component_left,component_right)
        
    def admittance(self):
        return 1/Add(self.component_left.admittance(),
            self.component_right.admittance())



class Component(Circuit):
    """docstring for Component"""
    def __init__(self, label):
        super(Component, self).__init__()
        self.label = label
    def remove_resistances(self):
        return self

class L(Component):
    def __init__(self, label):
        super(L,self).__init__(label)
    def admittance(self):
        return -sp.I*Mul(1/sp.Symbol('w'),1/sp.Symbol(self.label,real=True,nonnegative = True))

class J(L):
    def __init__(self, label):
        super(L,self).__init__(label)
        
class R(Component):
    def __init__(self, label):
        super(R, self).__init__(label)
    def admittance(self):
        return 1/sp.Symbol(self.label,real=True,nonnegative = True)

class C(Component):
    def __init__(self, label):
        super(C, self).__init__(label)
    def admittance(self):
        return sp.I*Mul(sp.Symbol('w'),sp.Symbol(self.label,real=True,nonnegative = True))
