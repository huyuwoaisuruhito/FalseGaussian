from timeit import default_timer as timer

import numpy as np
from scipy import linalg


def CI_A(HFArchive, level):
    CI( HFArchive.N, HFArchive.K,
        HFArchive.Hc, HFArchive.G, HFArchive.F,
        HFArchive.C, HFArchive.Vnn, 
        level, HFArchive)

def CI(N, K, Hc, G, F, C, Vnn, level, HFArchive=None):
    t = timer()
    print('\n========== Begin CI ==========')
    CI_C = []
    exc_list = []
    if level == 0 or level>min(N,K-N):
        __make_ci_c(C, N, K, CI_C, exc_list, min(N,K-N))
    else:
        __make_ci_c(C, N, K, CI_C, exc_list, level)
    CI_K = len(CI_C)
    t, dt = timer(), timer() - t
    print('Make {} excited determinant wavefunctions in {:.4f} s'.format(CI_K, dt))
 
    # n = 0
    # from multiprocessing import Pool
    # p = Pool(8)
    # def __callback(res):
    #     nonlocal CI_M, n, t
    #     if res != None:
    #         i, j = res[0]
    #         CI_M[i,j] = res[1]
    #         CI_M[j,i] = CI_M[i,j]
    #     n += 1
    #     if n % 1000 == 1:
    #         t, dt = timer(), timer() - t
    #         print('Computed '+ str(n-1) + ' of ' + str(CI_K) + '^2='
    #              + str(CI_K**2) + ' CI matrix unit. This cycle takes {:.4f} s.'.format(dt))

    CI_M = np.zeros((CI_K,CI_K))
    for i in range(CI_K):
        for j in range(CI_K):
                # p.apply_async(  __compute_ci_matrix_unit, 
                #                 args=(i, j, CI_C, N, K, Hc, G, CI_M),
                #                 callback=__callback)
                CI_M[i,j] = __compute_ci_matrix_unit(i, j, CI_C, N, K, Hc, G, F, CI_M)
    # p.close()
    # p.join()
    t, dt = timer(), timer() - t
    print('Compute {}x{} CI matrix in {:.4f} s'.format(CI_K, CI_K, dt))
    print(CI_M)
    print()
    if HFArchive:
        HFArchive.CI_M[level] = CI_M

    E = linalg.eigvalsh(CI_M)
    print(E)
    print('\nSolved CI matrix in {:.4f} s'.format(dt))
    print('====== CI in level:{} ======'.format(level))
    print('\nE = Eel + Vnn')
    print('  = {:.6f} + {:.6f}'.format(E[0], Vnn))
    print('  = {:.6f} Hartrees\n'.format(E[0]+Vnn))
    return E[0]

def __make_ci_c(C, N, K, CI_C, exc_list, max_c, count=0, m=0, n=0):
    CI_C.append(C)
    if count == max_c:
        return
    else:
        for i in range(m, N):
            for j in range(N+n, K):
                cC = np.copy(C)
                cC[:,[j,i]] = cC[:,[i,j]]
                __make_ci_c(cC, N, K, CI_C, max_c, count+1, i+1, j-N+1)

def __compute_ci_matrix_unit(n, m, CI_C, N, K, Hc, G, F, CI_M):
    # if abs(CI_M[n,m]) > 1e-10:
    #     return None
    C1, C2 = CI_C[n], CI_C[m]
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
                    F[i,j] += P[k,l] * (G[i,j,k,l] - 1/2 * G[i,k,j,l])
    
    E = 0.0
    for i in range(K):
        for j in range(K):
            E +=  1/2 * P[j,i] * (Hc[i,j] + F[i,j])
    print(F)
    return E