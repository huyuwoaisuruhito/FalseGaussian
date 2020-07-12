from timeit import default_timer as timer

import numpy as np
from scipy import linalg

def CI_A(HFArchive, level):
    CI( HFArchive.N, HFArchive.K,
        HFArchive.Hc, HFArchive.G, 
        HFArchive.C, HFArchive.Vnn, 
        level, HFArchive.hf_type, HFArchive.debug,HFArchive)

def CI(N, K, Hc, G, C, Vnn, level, hf_type=2, debug=0, HFArchive=None):
    t = timer()
    # N = N * hf_type
    # K = K * hf_type
    print('\n========== Begin CI ==========')

    exc_ob = [[[],[]]]
    if level == 0 or level>min(N,K-N):
        max_c = min(N,K-N)
    else:
        max_c = level
    __make_ci_c(N, K, max_c, n=N, res=exc_ob)
    exc_ob.sort(key=lambda x: len(x[0]))
    CI_K = len(exc_ob)

    t, dt = timer(), timer() - t
    print('Make {} excited determinant wavefunctions in {:.4f} s'.format(CI_K, dt))

    def __callback(res):
        nonlocal CI_M, n, t
        i, j = res[1:]
        CI_M[i,j] = res[0]
        CI_M[j,i] = CI_M[i,j]
        n += 1
        if n % 1000 == 1:
            t, dt = timer(), timer() - t
            print('Computed '+ str(n-1) + ' of ' + str(CI_K) + '^2='
                 + str(CI_K**2) + ' CI matrix unit. This cycle takes {:.4f} s.'.format(dt))
        if i!=j: n+=1

    n = 0
    CI_M = np.zeros((CI_K,CI_K))
    if CI_K > 30 and not debug:
        from multiprocessing import Pool
        p = Pool(8)
        for i, exc_ob1 in enumerate(exc_ob):
            for j, exc_ob2 in enumerate(exc_ob):
                if j <= i:
                    p.apply_async(  __compute_ci_matrix_unit, 
                                    args=(exc_ob1, exc_ob2, C, N, K, Hc, G, hf_type, i, j),
                                    callback=__callback)
        p.close()
        p.join()
    else:
        for i, exc_ob1 in enumerate(exc_ob):
            for j, exc_ob2 in enumerate(exc_ob):
                if j <= i:
                    CI_M[i,j] = __compute_ci_matrix_unit(exc_ob1, exc_ob2, C, N, K, Hc, G, hf_type, i, j)[0]
                    CI_M[j,i] = CI_M[i,j]

    t, dt = timer(), timer() - t
    print('Compute {}x{} CI matrix in {:.4f} s'.format(CI_K, CI_K, dt))
    print(CI_M)
    print()

    if HFArchive:
        HFArchive.CI_M[level] = CI_M

    e = linalg.eigvalsh(CI_M)

    print('e: ', e)
    print('Solved CI matrix in {:.4f} s'.format(dt))
    print('\n====== CI in level:{} ======'.format(level))
    print('\nE = Eel + Vnn')
    print('  = {:.6f} + {:.6f}'.format(e[0], Vnn))
    print('  = {:.6f} Hartrees\n'.format(e[0]+Vnn))
    return e[0]

def __make_ci_c(N, K, max_c, count=0, m=0, n=0, a=[], r=[], res=[]):
    if count==max_c:
        return
    for i in range(m, N):
        for j in range(n, K):
            a.append(i)
            r.append(j)
            res.append([a.copy(), r.copy()])
            __make_ci_c(N, K, max_c, count+1, i+1, j+1, a, r, res)
            a.pop()
            r.pop()

def __compute_ci_matrix_unit(exc_ob1, exc_ob2, C, N, K, Hc, G, hf_type, _i, _j):

    pC1 = np.array([i for i in range(K)])
    for a, r in zip(*exc_ob1):
        pC1[[r,a]] = pC1[[a,r]]
    
    pC2 = np.array([i for i in range(K)])
    for a, r in zip(*exc_ob2):
        pC2[[r,a]] = pC2[[a,r]]

    dC = [ (c1,c2) for c1,c2 in zip(pC1[:N],pC2[:N]) if c1!=c2 ] * hf_type
    
    unit = 0
    if len(dC) > 2:
        return unit, _i, _j
    
    if len(dC) == 2:
        m, p = dC[0]
        n, q = dC[0]
        _J, _K =0, 0
        for i in range(K):
            for j in range(K):
                for k in range(K):
                    for l in range(K):
                        if hf_type==2:
                            _J += C[i,m] * C[j,p] * C[k,n] * C[l,q] * G[i,j,k,l]
                        # _K -= C[i,m] * C[j,p] * C[k,n] * C[l,q] * G[i,k,j,l] / 2
        unit = _J + _K
        return unit, _i, _j
    
    C1 = np.copy(C)
    for a, r in zip(*exc_ob1):
        C1[:,[r,a]] = C1[:,[a,r]]

    C2 = np.copy(C)
    for a, r in zip(*exc_ob2):
        C2[:,[r,a]] = C2[:,[a,r]]

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
                    F[i,j] +=  P[k,l] * (G[i,j,k,l] - 1/2 * G[i,k,j,l])
    
    if len(dC) == 1:
        mu, nu = dC[0]
        for i in range(K):
            for j in range(K):
                unit += C[i,mu] * C[j,nu] * F[i,j]
        return unit, _i, _j

    if len(dC) == 0:
        for i in range(K):
            for j in range(K):
                unit +=  1/2 * P[j,i] * (Hc[i,j] + F[i,j])
        return unit, _i, _j