from timeit import default_timer as timer

import numpy as np
from scipy import linalg


def CI_A(HFArchive, level):
    CI( HFArchive.N, HFArchive.K,
        HFArchive.Hc, HFArchive.G, 
        HFArchive.C, HFArchive.Vnn, level)

def CI(N, K, Hc, G, C, Vnn, level):
    print(C)
    print()
    t = timer()
    print('\n========== Begin CI ==========')
    CI_C = []
    if level == 0:
        __make_ci_c(C, N, K, min(N,K-N), CI_C)
    else:
        __make_ci_c(C, N, K, level, CI_C)
    CI_K = len(CI_C)
    t, dt = timer(), timer() - t
    print('Make {} excited determinant wavefunctions in {:.4f} s'.format(CI_K, dt))

    CI_M = np.zeros((CI_K,CI_K))
    for i in range(CI_K):
        for j in range(CI_K):
            CI_M[i,j] = __compute_ci_matrix_unit(CI_C[i], CI_C[j], N, K, Hc, G)
    print('Compute {}x{} CI matrix in {:.4f} s'.format(CI_K, CI_K, dt))
    print(CI_M)

    E = linalg.eigvalsh(CI_M)
    print(E)
    print('\nSolved CI matrix in {:.4f} s'.format(dt))
    print('====== CI in level:{} ======'.format(level))
    print('\nE = Eel + Vnn')
    print('  = {:.6f} + {:.6f}'.format(E[0], Vnn))
    print('  = {:.6f} Hartrees\n'.format(E[0]+Vnn))
    return E[0]

def __make_ci_c(C, N, K, count, CI_C):
    CI_C.append(C)
    if count==0:
        return
    else:
        for i in range(N-count, N):
            for j in range(N, K):
                cC = np.copy(C)
                cC[:,[j,i]] = cC[:,[i,j]]
                __make_ci_c(cC, N, K, count-1, CI_C)

def __compute_ci_matrix_unit(C1, C2, N, K, Hc, G):
    P = np.zeros((K,K))
    for i in range(K):
        for j in range(K):
            for a in range(int(N)):
                P[i,j] += 2 * (C1[i,a] * C2[j,a])
    
    F = np.copy(Hc)
    for i in range(K):
        for j in range(K):
            for k in range(K):
                for l in range(K):
                    F[i,j] += P[k,l] * (G[i,j,k,l] - 1/2 * G[i,k,l,j])
    E = 0.0
    for i in range(K):
        for j in range(K):
            E +=  1/2 * P[j,i] * (Hc[i,j] + F[i,j])
    
    return E