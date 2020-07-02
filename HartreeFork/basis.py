import json
import numpy as np
import cmath

import HartreeFork.baseFunction as bf

root = 'HartreeFork/Basis/'

class Basis:
    def __init__(self, f):
        f_json = open(root + f, 'r').readlines()
        s_json = ''
        for s in f_json:
            s_json += s
        self.base = json.loads(s_json)
        self.basis = []

    def __get_base(self, atom_number, R):
        basis = []
        ls = self.base["elements"][str(atom_number)]["electron_shells"]
        for l in ls:
            if l["function_type"] == "gto":
                for i, ang in enumerate(l["angular_momentum"]):
                    def _(ang, l_ang, i):
                        if ang == 0:
                            _len = len(l["exponents"])
                            cgf = bf.CGF(k=1, R=R, ang=l_ang, *zip([float(k) for k in l["coefficients"][i]], [float(exp) for exp in l["exponents"]]))
                            basis.append(cgf)
                        else:
                            for j in range(3):
                                l_ang[j] += 1
                                _(ang-1, l_ang, i)
                                l_ang[j] -= 1
                    _(ang, [0, 0, 0], i)
        return basis

    def S(self):
        _len = len(self.basis)
        S = np.zeros((_len,_len))
        for i in range(_len):
            for j in range(i, _len):
                s = (self.basis[i] * self.basis[j]).CGF_integration()
                S[i][j] = s
                S[j][i] = s
        return np.matrix(S)


    def __symmetric_orthogonalization(self):
        S = self.S()
        print(S)
        s, U = np.linalg.eigh(S)
        s = np.diag([1/cmath.sqrt(ss) for ss in s])
        X = np.matmul(U, s)
        return X

    def get_basis(self, atoms):
        print('#CreatBasis')
        for n, x, y, z in atoms:
            self.basis += self.__get_base(n, np.array([x, y, z]))
        print('#SymmetricOrthogonalizations')
        X = self.__symmetric_orthogonalization()
        return (X, self.basis)


class Basis_3_21g (Basis):
    def __init__(self):
        super().__init__('3-21g.json')

class Basis_STO_3g (Basis):
    def __init__(self):
        super().__init__('sto-3g.json')



atoms = [[6, -1.2e-07, 0.27586206, 0.0], [1, -1.07000012, 0.27586206, 0.0], [7, 1.14659988, 0.27586206, 0.0]]
basis = Basis_3_21g()
basis.get_basis(atoms)
print('Clear')