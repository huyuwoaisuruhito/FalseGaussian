from timeit import default_timer as timer

import numpy as np
from scipy import linalg

Rate = 0.6

def HF_A(HFArchive):
    return HF(  HFArchive.N, HFArchive.K, 
                HFArchive.S, HFArchive.Hc, 
                HFArchive.G, HFArchive.Vnn, 
                HFArchive.hf_type, HFArchive.SCF_MAX_iteration, 
                HFArchive.SCF_ERROR, HFArchive.debug, 
                HFArchive )

def HF(N, K, S, Hc, G, Vnn, oS, SCF_MAX_iteration=200, SCF_ERROR=1e-6, debug=0, HFArchive=None):
    '''
    oS = 1 RHF
    oS = 2 UHF
    '''
    print('\n========== Begin {} =========='.format(['RHF', 'UHF'][oS-1]))
    
    if oS == 2 and len(N)==1:
        N = [ N[0]-N[0]//2, N[0]//2 ]
    elif oS == 1 and len(N)==1:
        N = [N[0]//2]
    else:
        raise RuntimeError()

    print('\nElectrons: '+str(N)+'\n')

    X = linalg.sqrtm(linalg.inv(S))
    C = np.zeros((oS,K,K))
    P = np.zeros((oS,K,K))

    inisial_gasse(oS, K, P)

    E = 0.0
    count = 0

    t = timer()
    for iteration in range(SCF_MAX_iteration):
        E_old = E
        E = 0.0
        F = np.zeros((oS,K,K))
        nP = np.zeros((oS,K,K))
        for s in range(oS):
            nE = 0.0
            for i in range(K):
                for j in range(K):
                    F[s,i,j] = Hc[i,j]
                    for k in range(K):
                        for l in range(K):
                            F[s,i,j] -=  1/2  * P[s,k,l] * G[i,k,j,l]
                            for t in range(oS):
                                F[s,i,j] +=  1/oS * P[t,k,l] * G[i,j,k,l]

            Fp = X.T @ F[s] @ X
            e, Cp = linalg.eigh(Fp)
            C[s] = X @ Cp

            for i in range(K):
                for j in range(K):
                    for a in range(N[s]):
                        nP[s,i,j] += 2 * (C[s,i,a] * C[s,j,a])
            
            P[s] = nP[s]*Rate + P[s]*(1-Rate)

            if P.any()>2:
                print('Warning: Density Matrix Overflow')
            # print('\nDensity Matrix (iteration {:2d})'.format(count))
            # print(str(P))

            for i in range(K):
                for j in range(K):
                    nE +=  1/2 * P[s,j,i] * (Hc[i,j] + F[s,i,j]) / oS
            # print('s: ', s, 'E: ', nE)
            E += nE
         
        t, dt = timer(), timer() - t
        print('E (iteration {:2d}) = {:12.6f} \t in {:.4f} s'.format(count,E, dt))
        if (abs(E-E_old) < SCF_ERROR) and (iteration>0):
            print('\n====== SCP converged in {} steps ======'.format(count))
            print('\nE = Eel + Vnn')
            print('  = {:.6f} + {:.6f}'.format(E, Vnn))
            print('  = {:.6f} Hartrees ({} iterations)\n'.format(E+Vnn, count))
            # P_ = np.zeros((K,K))
            # for s in range(oS):
            #     P_ += P[s,:,:]
            # print( P_ )
            # print('e: ', e)
            # print(str(P))
            print()
            if HFArchive:
                HFArchive.e = e
                HFArchive.X = X
                HFArchive.C = C
                HFArchive.P = P
            return E + Vnn
        
        count += 1
    
    print('!!!SCF iteration does not converge!!!')
    return -1

def inisial_gasse(oS, K, P):
    for s in range(oS):
        np.random.seed(10)
        P[s,:,:] = np.random.rand(K,K)