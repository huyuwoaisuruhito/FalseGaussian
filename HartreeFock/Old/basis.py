import json
import numpy as np
import cmath

import baseFunction as bf

root = 'HartreeFork/Basis/'

class _Basis:
    def __init__(self, file, flag=1):

        #TODO:slater orbital exponents;flag

        f_json = open(root + file, 'r').readlines()
        s_json = ''
        for s in f_json:
            s_json += s
        self.base = json.loads(s_json)["elements"]
        self.basis = []
        self.flag = flag

    def __make_base(self, atom_number, R):
        basis = []
        ls = self.base[str(atom_number)]["electron_shells"]
        for l in ls:
            if l["function_type"] == "gto":
                for i, ang in enumerate(l["angular_momentum"]):
                    def _(ang, l_ang, i):
                        if ang == 0:
                            _len = len(l["exponents"])
                            if self.flag:
                                cgf = bf.CGF(k=1, R=R, ang=l_ang, *zip([float(k) for k in l["coefficients"][i]], [float(exp) for exp in l["exponents"]]))
                            else:
                                cgf = bf.CGF(k=1, R=R, ang=l_ang, *zip([float(k) for k in l["coefficients"][i]], [float(exp)*atom_number**2 for exp in l["exponents"]]))
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
                s = (self.basis[i] * self.basis[j]).CGF_integration(1)
                S[i][j] = s
                S[j][i] = s
        print('S:')
        print(S)
        return np.matrix(S)

    def __symmetric_orthogonalization(self):
        S = self.S()
        s, U = np.linalg.eigh(S)
        s = np.diag([1/cmath.sqrt(ss) for ss in s])
        X = np.matmul(U, s)
        return X

    def make_basis(self, atoms, type='ortho'):
        X = []
        print('#CreatBasis')
        for n, x, y, z in atoms:
            self.basis += self.__make_base(n, np.array([x, y, z]))
        print('#SymmetricOrthogonalizations')
        if type == 'ortho':
            X = self.__symmetric_orthogonalization()
        return (self.basis, X)


class Basis_3_21g(_Basis):
    def __init__(self):
        super().__init__('3-21g.json')

class Basis_STO_3g(_Basis):
    def __init__(self):
        super().__init__('sto-3g.json')

class Basis_STO_3g_test(_Basis):
    def __init__(self):
        super().__init__('sto_3g_test.json')