from timeit import default_timer as timer

import numpy as np
from scipy import linalg

Rate = 0.6

def HF_A(HFArchive):
    RHF(HFArchive.N, HFArchive.K, 
        HFArchive.S, HFArchive.Hc, 
        HFArchive.G, HFArchive.Vnn, 
        HFArchive.SCF_MAX_iteration, 
        HFArchive.SCF_ERROR, HFArchive)

def RHF(N, K, S, Hc, G, Vnn, SCF_MAX_iteration=200, SCF_ERROR=1e-6, HFArchive=None):
    X = linalg.sqrtm(linalg.inv(S))
    C = np.zeros((K,K))
    P = np.zeros((K,K))
    E = 0.0
    count = 0

    print('\n====== Begin the SCF iteration ======')
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
                for a in range(int(N[0])):
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
        if (abs(E-E_old) < SCF_ERROR) and (iteration>0):
            print('\n====== SCP converged in {} steps ======'.format(count))
            print('\nE = Eel + Vnn')
            print('  = {:.6f} + {:.6f}'.format(E, Vnn))
            print('  = {:.6f} Hartrees ({} iterations)\n'.format(E+Vnn, count))
            print('e: ', e)
            print()
            if HFArchive:
                HFArchive.e = e
                HFArchive.X = X
                HFArchive.C = C
                HFArchive.P = P
                HFArchive.F = F
            return E + Vnn
        
        count += 1
    
    print('!!!SCF iteration does not converge!!!')
    return -1

