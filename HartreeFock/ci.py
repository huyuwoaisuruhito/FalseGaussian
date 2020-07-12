from timeit import default_timer as timer

import numpy as np
from scipy import linalg

def CI_A(HFArchive, level):
    return CI(  HFArchive.N, HFArchive.K,
                HFArchive.Hc, HFArchive.G, 
                HFArchive.C, HFArchive.Vnn, 
                level, HFArchive.hf_type, 
                HFArchive.debug, HFArchive )

def CI(N, K, Hc, G, C, Vnn, level, oS=1, debug=0, HFArchive=None):
    '''
    CI方法，已测试
    level: 计算到第level重激发，为0时计算FCI
    '''
    t = timer()
    print('\n========== Begin CI ==========\n')
    
    if oS == 2 and len(N)==1:
        N = [ N[0]-N[0]//2, N[0]//2 ]
    elif oS == 1 and len(N)==1:
        N = [N[0]//2]
    else:
        raise RuntimeError()
    
    sN = sum(N)
    exc_ob = [[[],[]]]
    if level == 0 or level > min(sN, K-sN):
        max_c = min(sN, K-sN)
    else:
        max_c = level
    __make_ci_c(sN, K, max_c, n=sN, res=exc_ob)
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
                                    args=(exc_ob1, exc_ob2, C, N, K, Hc, G, oS, i, j),
                                    callback=__callback)
        p.close()
        p.join()
    else:
        for i, exc_ob1 in enumerate(exc_ob):
            for j, exc_ob2 in enumerate(exc_ob):
                if j <= i:
                    CI_M[i,j] = __compute_ci_matrix_unit(exc_ob1, exc_ob2, C, N, K, Hc, G, oS, i, j)[0]
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
    
    return e[0]+Vnn


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


def __compute_ci_matrix_unit(exc_ob1, exc_ob2, C, N, K, Hc, G, oS, _i, _j):
    '''
    oS = 1
    构造RHF的CI矩阵，已测试
    '''
    '''
    oS = 2
    构造UHF的CI矩阵，这是毫无必要的，，，虽然我写了，但这里还有问题。
    在UHF下，N = [Na, Nb]，K为总自旋轨道数目，重新定义K为总RHF轨道数目
    '''

    pC1 = np.array([i for i in range(K)])
    for a, r in zip(*exc_ob1):
        pC1[[r,a]] = pC1[[a,r]]
    
    pC2 = np.array([i for i in range(K)])
    for a, r in zip(*exc_ob2):
        pC2[[r,a]] = pC2[[a,r]]

    sN = sum(N)
    dC = [ (c1%oS, c1//oS, c2%oS, c2//oS) for c1,c2 in zip(pC1[:sN],pC2[:sN]) if c1!=c2 ] * (3 - oS)
    
    unit = 0
    if len(dC) > 2:
        return unit, _i, _j
    
    if len(dC) == 2:
        sm, m, sp, p = dC[0]
        sn, n, sq, q = dC[1]
        if not (sm == sp and sn == sq) and not (sm == sn and sp == sq):
            return unit, _i, _j
        for i in range(K):
            for j in range(K):
                for k in range(K):
                    for l in range(K):
                        if sm == sp and sn == sq:
                            unit += C[sm,i,m] * C[sp,j,p] * C[sn,k,n] * C[sq,l,q] * G[i,j,k,l] #[mp|nq]=J
                        if sm == sn and sp == sq and oS==2:
                            unit += C[sm,i,m] * C[sn,j,n] * C[sp,k,p] * C[sq,l,q] * G[i,j,k,l] #[mn|pq]=K
        return unit, _i, _j

    P = np.zeros((oS,K,K))
    for s in range(oS):

        C1 = np.copy(C)
        for a, r in zip(*exc_ob1):
            C1[s,:,[r,a]] = C1[s,:,[a,r]]

        C2 = np.copy(C)
        for a, r in zip(*exc_ob2):
            C2[s,:,[r,a]] = C2[s,:,[a,r]]
        
        for i in range(K):
            for j in range(K):
                for a in range(N[s]):
                    P[s,i,j] += 2 * (C1[s,i,a] * C2[s,j,a])

    F = np.zeros((oS,K,K))
    for s in range(oS):
        for i in range(K):
            for j in range(K):
                F[s,i,j] = Hc[i,j]
                for k in range(K):
                    for l in range(K):
                        F[s,i,j] -= 1/2 * P[s,k,l] * G[i,k,j,l]
                        for t in range(oS):
                            F[s,i,j] += 1/oS * P[t,k,l] * G[i,j,k,l]
    
    if len(dC) == 1:
        sm, m, sp, p = dC[0]
        if sm != sp:
            return unit, _i, _j
        for i in range(K):
            for j in range(K):
                unit += C1[sm,i,m] * C2[sp,j,p] * F[s,i,j]
        return unit, _i, _j

    if len(dC) == 0:
        for s in range(oS):
            for i in range(K):
                for j in range(K):
                    unit +=  1/2 * P[s,j,i] * (Hc[i,j] + F[s,i,j]) / oS
        return unit, _i, _j