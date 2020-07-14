import os.path
import pickle
from timeit import default_timer as timer

import numpy as np

import HartreeFock.basis as bs
import HartreeFock.ci as ci
import HartreeFock.molecularIntegrals as mi
import HartreeFock.hf as hf

np.set_printoptions(threshold=np.inf)
np.set_printoptions(precision=5)
np.set_printoptions(linewidth=2000)
np.set_printoptions(suppress=True)

A0 = 0.52917721067


class HFArchive:

    def __init__(self, N, atoms, bname, fname, hf_type=1, expand=0, fresh=0):

        self.TEMP_DIR = '.temp'
        self.SCF_MAX_iteration = 1000
        self.SCF_ERROR = 1e-8
        self.debug = 0

        self.N = [N]
        A_2_atom_unit(atoms)
        self.atoms = atoms
        self.expand = expand
        self.fresh = fresh
        path = '{0}/{1}'.format(self.TEMP_DIR, fname)
        self.name = path + '/{0}.{1}'.format(fname, bname)
        if not os.path.exists(path):
            os.makedirs(path)

        self.basis = bs.Basis(bname+'.json').make_basis(atoms, expand)
        self.K = len(self.basis)

        self.Vnn = None
        self.S = None
        self.T = None
        self.V = None
        self.G = None

        self.X = None
        self.C = None
        self.P = None
        self.H = None

    def init_molecular_integrals(self):
        self.Vnn = mi.nuclear_repulsion(self.atoms)
        self.S, self.T, self.V, self.G = make_molecular_integrals(
            self.K, self.basis, self.atoms, self.name, self.fresh)
        self.Hc = self.T + self.V

    def dump(self):
        with open(self.name+'.arc', 'wb') as file:
            pickle.dump(self, file)

    @staticmethod
    def load(fname):
        with open(fname+'.arc', 'rb') as file:
            return pickle.load(file)

    def RHF(self):
        return hf.HF_A(self, 1)

    def UHF(self):
        return hf.HF_A(self, 2)

    def CI(self, level):
        return ci.CI_A(self, level)

    def FCI(self):
        return ci.CI_A(self, 0)


def A_2_atom_unit(atoms):
    for a in atoms:
        for i in range(1, 4):
            a[i] /= A0


def make_molecular_integrals(K, basis, atoms, name, fresh):
    t = timer()
    if os.path.exists(name+'.S.npy') and not fresh:
        S = np.load(name+'.S.npy')
    else:
        S = np.zeros((K, K))
        mi.buildS(basis, S)
        dump_matrix(name, 'S', S)
    t, dt = timer(), timer() - t
    print('Prepare overlap matrix in:{:.4f} s'.format(dt))

    if os.path.exists(name+'.T.npy') and not fresh:
        T = np.load(name+'.T.npy')
    else:
        T = np.zeros((K, K))
        mi.buildT(basis, T)
        dump_matrix(name, 'T', T)
    t, dt = timer(), timer() - t
    print('Prepare kinetic matrix in:{:.4f} s'.format(dt))

    if os.path.exists(name+'.V.npy') and not fresh:
        V = np.load(name+'.V.npy')
    else:
        V = np.zeros((K, K))
        mi.buildV(basis, atoms, V)
        dump_matrix(name, 'V', V)
    t, dt = timer(), timer() - t
    print('Prepare nuclear attraction matrix in:{:.4f} s'.format(dt))

    if os.path.exists(name+'.G.npy') and not fresh:
        G = np.load(name+'.G.npy')
    else:
        G = np.zeros((K, K, K, K))
        if K > 8:
            mi.buildG_P(basis, G)
        else:
            mi.buildG(basis, G)
        dump_matrix(name, 'G', G)
    t, dt = timer(), timer() - t
    print('Prepare two-electron integral matrix in:{:.4f} s'.format(dt))

    return (S, T, V, G)


def dump_matrix(name, N, M, txt=1):
    np.save(open('{}.{}.npy'.format(name, N), 'wb'), M)
    if txt:
        open('{}.{}.txt'.format(name, N), 'w').write(str(M))
