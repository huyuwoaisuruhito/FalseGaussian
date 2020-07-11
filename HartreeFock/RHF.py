import os.path
from timeit import default_timer as timer

import numpy as np
from scipy import linalg

import HartreeFock.basis2 as bs2
import HartreeFock.molecularIntegrals as mi

debug = 0

SCF_MAX_iteration = 200
SCF_ERROE = 1e-6
TEMP_DIR = '.temp'
Rate = 0.8

np.set_printoptions(threshold=np.inf)
np.set_printoptions(precision=5)
np.set_printoptions(linewidth=200)
np.set_printoptions(suppress=True)

def RHF(N, atoms, bname, fname):
    name = '{}/{}.{}'.format(TEMP_DIR, fname, bname)

    B = bs2.Basis(bname+'.json')
    basis = B.make_basis(atoms)

    K = len(basis)
    S, Hc, G = make_molecular_integrals(K, basis, atoms, name)
    Vnn = mi.nuclear_repulsion(atoms)

    X = linalg.sqrtm(linalg.inv(S))
    C = np.zeros((K,K))
    P = np.zeros((K,K))
    E = 0.0
    count = 0

    print('\nBegin the SCF iteration:')
    t = timer()
    for iteration in range(SCF_MAX_iteration):
        E_old = E
        F = np.copy(Hc)
        for i in range(K):
            for j in range(K):
                for k in range(K):
                    for l in range(K):
                        F[i,j] += P[k,l] * (G[i,j,k,l] - 1/2 * G[i,k,l,j])

        Fp = X.T @ F @ X
        e, Cp = linalg.eigh(Fp)
        C = X @ Cp

        nP = np.zeros((K,K))
        for i in range(K):
            for j in range(K):
                for a in range(int(N/2)):
                    nP[i,j] += 2 * (C[i,a] * C[j,a])
        
        if (iteration>0):
            P = nP*Rate + P*(1-Rate)
        
        if P.any()>2:
            print('Warning: Density Matrix Overflow')
        # print('Density Matrix (iteration {:2d})'.format(count))
        # print('\n'+str(P))

        E = 0.0
        for i in range(K):
            for j in range(K):
                E += 1/2 * P[j,i] * (Hc[i,j] + F[i,j])
        
        t, dt = timer(), timer() - t
        print('E (iteration {:2d}) = {:12.6f} \t in {:.4f} s'.format(count,E, dt))
        if (abs(E-E_old) < SCF_ERROE) and (iteration>0):
            print('\n====== SCP converged in {} steps ======'.format(count))
            print('\nE = Eel + Vnn')
            print('  = {:.6f} + {:.6f}'.format(E, Vnn))
            print('  = {:.6f} Hartrees ({} iterations)\n'.format(E+Vnn, count))
            dump_matrix(name, 'C', C)
            dump_matrix(name, 'P', P)
            return E + Vnn
        
        count += 1
    
    print('!!!Iteration does not converge!!!')
    return -1

def make_molecular_integrals(K, basis, atoms, name):
    
    t = timer()
    if os.path.exists(name+'.S.npy') and not debug:
        S = np.load(name+'.S.npy')
    else:
        S = np.zeros((K,K))
        mi.buildS(basis, S)
        dump_matrix(name, 'S', S)
    t, dt = timer(), timer() - t
    print('Prepare overlap matrix in:{:.4f} s'.format(dt))

    if os.path.exists(name+'.T.npy') and not debug:
        T = np.load(name+'.T.npy')
    else:
        T = np.zeros((K,K))
        mi.buildT(basis, T)
        dump_matrix(name, 'T', T)
    t, dt = timer(), timer() - t
    print('Prepare kinetic matrix in:{:.4f} s'.format(dt))

    if os.path.exists(name+'.V.npy') and not debug:
        V = np.load(name+'.V.npy')
    else:
        V = np.zeros((K,K))
        mi.buildV(basis, atoms, V)
        dump_matrix(name, 'V', V)
    t, dt = timer(), timer() - t
    print('Prepare nuclear attraction matrix in:{:.4f} s'.format(dt))

    if os.path.exists(name+'.G.npy') and not debug:
        G = np.load(name+'.G.npy')
    else:
        G = np.zeros((K,K,K,K))
        mi.buildG_P(basis, G)
        dump_matrix(name, 'G', G)
    t, dt = timer(), timer() - t
    print('Prepare two-electron integral matrix in:{:.4f} s'.format(dt))

    Hc = T + V
    return (S, Hc, G)

def dump_matrix(name, N, M):
    np.save(open('{}.{}.npy'.format(name, N),'wb'), M)
    open('{}.{}.txt'.format(name, N),'w').write( str(M) )

if __name__ == '__main__':
    atoms = [[8, 0.000000, 0.000000, 0.227000], [1, 0.000000, 1.353000,-0.908000], [1, 0.000000,-1.353000,-0.908000]]
    RHF(10, atoms, '3-21g', 'H2O')
