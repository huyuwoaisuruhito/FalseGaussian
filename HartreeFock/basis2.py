import json
import numpy as np
from scipy import misc, special

root = 'HartreeFock/Basis/'

class PGF:
    def __init__(self, k, exp, R, ang):
        self.k = k
        self.exp = exp
        self.R = np.array(R)
        self.ang = ang

class CGF:
    def __init__(self, k, R, ang, *arg):
        self.k = k
        self.cgf = []
        for a in arg:
            self.cgf.append(PGF(*a, R, ang))

class Basis:
    def __init__(self, file):
        f_json = open(root + file, 'r').readlines()
        s_json = ''
        for s in f_json:
            s_json += s
        self.base = json.loads(s_json)["elements"]
        self.basis = []

    def __make_base(self, atom_number, R, expand=0):
        basis = []
        ls = self.base[str(atom_number)]["electron_shells"]
        for l in ls:
            if l["function_type"] == "gto":
                for i, ang in enumerate(l["angular_momentum"]):
                    def _(ang, l_ang, i):
                        if ang == 0:
                            co = l["coefficients"][i]
                            if expand:
                                for k, exp in zip([float(k) for k in co], [float(exp) for exp in l["exponents"]]):
                                    cgf = CGF(1, R, l_ang.copy(), (1, exp))
                                    basis.append(cgf)
                            else:
                                cgf = CGF(1, R, l_ang.copy(), *zip([float(k) for k in co], [float(exp) for exp in l["exponents"]]))
                                basis.append(cgf)
                        else:
                            for j in range(3):
                                l_ang[j] += 1
                                _(ang-1, l_ang, i)
                                l_ang[j] -= 1
                    _(ang, [0, 0, 0], i)
        return basis

    def make_basis(self, atoms, type='ortho', expand=0):
        for n, x, y, z in atoms:
            self.basis += self.__make_base(n, np.array([x, y, z]), expand)
        return self.basis

class Basis_STO_3g(Basis):
    def __init__(self):
        super().__init__('sto-3g.json')

class Basis_3_21g(Basis):
    def __init__(self):
        super().__init__('3-21g.json')

class Basis_6_31g(Basis):
    def __init__(self):
        super().__init__('6-31g.json')