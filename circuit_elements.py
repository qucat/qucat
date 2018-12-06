import lcapy

class L(lcapy.oneport.L):
    def __init__(self, label):
        super(L,self).__init__(label)
        self.label = label

class J(lcapy.oneport.L):
    def __init__(self, label):
        super(J,self).__init__(label)
        self.label = label

class R(lcapy.oneport.R):
    def __init__(self, label):
        super(R, self).__init__(label)
        self.label = label

class C(lcapy.oneport.C):
    def __init__(self, label):
        super(C, self).__init__(label)
        self.label = label