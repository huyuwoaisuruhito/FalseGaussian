import numpy as np

import basis2 as bs2


class InterfaceHF:
    def __init__(self, mainwindow, molecule):
        self._parent = mainwindow
        self._molecule = molecule


class ComputingRHF:
    def __init__(self, atoms, object_basis, electrons):
        self.electrons = electrons
        self.object_basis = object_basis
        self.atoms = atoms
        self.B = None
        self.X = None
        self.P = None
        self.F = None
        self.XFX = None

    def run(self):
        self.geometry()
        self.make_basis()
        self.compute_integration()

    def geometry(self):
        pass

    def make_basis(self):
        self.B, self.X = self.object_basis.make_basis(self.atoms)

    def compute_integration(self):
        self.__O1()

    def __O1(self):
        n = len(self.B)

        T = np.zeros((n, n))
        for i in range(n):
            phi_1 = self.B[i]
            for j in range(n):
                phi_2 = self.B[j]
                d_phi_2 = phi_2.CGF_diff()
                for k in d_phi_2:
                    T[i][j] += (phi_1 * k).CGF_integration()
                T[i][j] *= -1/2
        print('T:\n', T)

        V = np.zeros((n, n))
        for Z, ax, ay, az in self.atoms:
            for i in range(n):
                phi_1 = self.B[i]
                for j in range(n):
                    phi_2 = self.B[j]
                    V[i][j] += -Z
        return T + V

    def __O2(self):
        pass

    def initial_guess_of_P(self):
        pass



atoms = [[1, 0.0, 0.0, 0.0], [2, 0.0, 0.0, 1.4632]]
# object_basis = bs.Basis_STO_3g_test()
object_basis = bs.Basis_STO_3g()
RHF = ComputingRHF(atoms, object_basis, 2)
RHF.run()
print('Clear')
