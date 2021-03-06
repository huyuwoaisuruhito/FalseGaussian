from timeit import default_timer as timer

import numpy as np
from scipy import linalg


def CI_A(HFArchive, level):
    return CI(HFArchive.N, HFArchive.K, HFArchive.X, 
              HFArchive.Hc, HFArchive.G,
              HFArchive.C, HFArchive.Vnn, HFArchive.F,
              level, 2,
              HFArchive.debug, HFArchive)


def CI(N, K, X, Hc, G, C, Vnn, _F, level, oS=1, debug=0, HFArchive=None):
    '''
    CI方法，已测试
    level: 计算到第level重激发，为0时计算FCI
    '''
    t = timer()
    print('\n========== Begin CI ==========\n')

    if oS == 1 and len(N) == 1:
        N = [N[0]//2]
    elif oS == 2:
        if len(N) == 1:
            N = [N[0]-N[0]//2, N[0]//2]
        elif len(N) == 2:
            pass
        else:
            raise RuntimeError()
    else:
        raise RuntimeError()

    exc_ob = [[]]*3
    for s in range(oS):
        exc_ob_s = [[[], []]]
        if level == 0 or level > min(N[s], K-N[s]):
            max_c = min(N[s], K-N[s])
        else:
            max_c = level
        __make_ci_c(N[s], s, oS, K, max_c, n=N[s], res=exc_ob_s)
        exc_ob[s] = exc_ob_s[1:]
        if s == 1:
            for n in exc_ob[0]:
                for m in exc_ob[1]:
                    zipped = [n[0]+m[0], n[1]+m[1]]
                    exc_ob[2].append(zipped)
    exc_ob = [[[], []]] + exc_ob[0] + exc_ob[1] + exc_ob[2]
    exc_ob.sort(key=lambda x: len(x[0]))
    if level: exc_ob = [ _ for _ in exc_ob if len(_[0])<=level]
    print(exc_ob)
    CI_K = len(exc_ob)

    t, dt = timer(), timer() - t
    print('Make {} excited determinant wavefunctions in {:.4f} s\n'.format(CI_K, dt))

    E = __naive_method(sum(N), N, K, X, Hc, G, C, Vnn, _F, 
                       level, oS, CI_K, exc_ob, debug)

    print('Solved CI matrix in {:.4f} s'.format(dt))
    print('\n====== CI in level:{} ======'.format(level))
    print('\nE = Eel + Vnn')
    print('  = {:.6f} + {:.6f}'.format(E, Vnn))
    print('  = {:.6f} Hartrees\n'.format(E+Vnn))

    return E + Vnn


def __naive_method(sN, N, K, X, Hc, G, C, Vnn, _F, level, oS, CI_K, exc_ob, debug):

    def __callback(res):
        nonlocal CI_M, n, t
        i, j = res[1:3]
        CI_M[i, j] = res[0]
        CI_M[j, i] = CI_M[i, j]
        n += 1
        if n % 1000 == 1:
            t, dt = timer(), timer() - t
            print('Computed ' + str(n-1) + ' of ' + str(CI_K) + '^2='
                  + str(CI_K**2) + ' CI matrix unit. This cycle takes {:.4f} s.'.format(dt))
        if i != j:
            n += 1

    n = 0
    t = timer()
    CI_M = np.zeros((CI_K, CI_K))
    if CI_K > 30 and not debug:
        from multiprocessing import Pool
        p = Pool(8)
        for i, exc_ob1 in enumerate(exc_ob):
            for j, exc_ob2 in enumerate(exc_ob):
                if j <= i:
                    p.apply_async(__compute_ci_matrix_unit,
                                  args=(exc_ob1, exc_ob2, C, sN,
                                        N, K, X, Hc, G, oS, i, j, _F),
                                  callback=__callback)
        p.close()
        p.join()
    else:
        for i, exc_ob1 in enumerate(exc_ob):
            for j, exc_ob2 in enumerate(exc_ob):
                if j <= i:
                    CI_M[i, j] = __compute_ci_matrix_unit(
                        exc_ob1, exc_ob2, C, sN, N, K, X, Hc, G, oS, i, j, _F)[0]
                    CI_M[j, i] = CI_M[i, j]

    t, dt = timer(), timer() - t
    print('Compute {}x{} CI matrix in {:.4f} s'.format(CI_K, CI_K, dt))
    print(CI_M)
    print()

    e ,C = linalg.eigh(CI_M)
    print('\ne: \n', e)
    print('\nC:\n ', C)

    return e[0]


def __make_ci_c(N, s, oS, K, max_c, count=0, m=0, n=0, a=[], r=[], res=[]):
    if count == max_c:
        return
    for i in range(m, N):
        for j in range(n, K):
            a.append(i*oS+s)
            r.append(j*oS+s)
            res.append([a.copy(), r.copy()])
            __make_ci_c(N, s, oS, K, max_c, count+1, i+1, j+1, a, r, res)
            a.pop()
            r.pop()


def __compute_ci_matrix_unit(exc_ob1, exc_ob2, C, sN, N, K, X, Hc, G, oS, _i, _j, _F):
    '''
    oS = 1
    构造RHF的CI矩阵，已测试
    '''
    '''
    oS = 2
    构造UHF的CI矩阵，这是毫无必要的，，，虽然我写了，但这里还有问题。
    在UHF下，N = [Na, Nb]，K为总自旋轨道数目，重新定义K为总RHF轨道数目
    '''
    if C.shape[0] == 1:
        C = np.append(C, np.copy(C), axis=0)
        _F = np.append(_F, np.copy(_F), axis=0)

    pC1 = np.array([i for i in range(K*oS)])
    for a, r in zip(*exc_ob1):
        a, r = a, r
        pC1[[r, a]] = pC1[[a, r]]

    pC2 = np.array([i for i in range(K*oS)])
    for a, r in zip(*exc_ob2):
        a, r = a, r
        pC2[[r, a]] = pC2[[a, r]]

    dC = [(c1//oS, c1 % oS, c2//oS, c2 % oS)
          for c1, c2 in zip(pC1[:sN], pC2[:sN]) if c1 != c2] * (3 - oS)

    active = pC1[:sN] + pC2[:sN]

    unit = 0
    if len(dC) > 2:# or ( (_i == 0 or _j == 0) and len(dC) == 1 ):
        return unit, _i, _j

    if len(dC) == 2:
        m, sm, p, sp = dC[0]
        n, sn, q, sq = dC[1]
        if not (sm == sp and sn == sq) and not (sm == sn and sp == sq):
            return unit, _i, _j
        for i in range(K):
            for j in range(K):
                for k in range(K):
                    for l in range(K):
                        if sm == sp and sn == sq:
                            unit += (C[sm, i, m] * C[sp, j, p] * C[sn, k, n]
                                     * C[sq, l, q] * G[i, j, k, l])  # [mp|nq]=J
                        if sm == sn and sp == sq and oS == 2:
                            unit -= (C[sm, i, m] * C[sn, j, n] * C[sp, k, p]
                                     * C[sq, l, q] * G[i, j, k, l])  # [mn|pq]=K
        return unit, _i, _j

    C1 = np.copy(C)
    C2 = np.copy(C)
    for a, r in zip(*exc_ob1):
        assert(a % oS == r % oS)
        sa, a, r = a % oS, a//oS, r//oS
        C1[sa, :, [r, a]] = C1[sa, :, [a, r]]

    for a, r in zip(*exc_ob2):
        assert(a % oS == r % oS)
        sa, a, r = a % oS, a//oS, r//oS
        C2[sa, :, [r, a]] = C2[sa, :, [a, r]]

    P = np.zeros((oS, K, K))
    for i in range(K):
        for j in range(K):
            for s in range(oS):
                for a in range(N[s]):
                    na1 = pC1[2*a+s]//2
                    na2 = pC2[2*a+s]//2
                    if na1 == na2:
                        P[s, i, j] += C[s, i, na1] * C[s, j,na2]

    F = np.zeros((oS, K, K))
    for s in range(oS):
        for i in range(K):
            for j in range(K):
                F[s, i, j] = Hc[i, j]
                for k in range(K):
                    for l in range(K):
                            F[s, i, j] -= P[s, k, l] * G[i, k, j, l]
                            F[s, i, j] += P[0, k, l] * G[i, j, k, l]
                            F[s, i, j] += P[1, k, l] * G[i, j, k, l]
    
    if len(dC) == 1:
        I = np.zeros((oS, K, K))
        for s in range(oS):
            for i in range(K):
                for j in range(K):
                    for k in range(K):
                        for l in range(K):
                            I[s, k, l] += C[s, i, k] * C[s, j, l] * F[s, i, j]
        # print(dC[0][1], dC[0][0], dC[0][2])
        # print(I[dC[0][1], dC[0][0], dC[0][2]])
        # print(I)
        # print()

    if len(dC) == 1:
        m, sm, p, sp = dC[0]
        if sm != sp:
            return unit, _i, _j
        for i in range(K):
            for j in range(K):
                unit += C[sm, i, m] * C[sp, j, p] * F[sm, i, j]
        return unit, _i, _j

    if len(dC) == 0:
        for s in range(oS):
            for i in range(K):
                for j in range(K):
                    unit += 1/2 * P[s, j, i] * (Hc[i, j] + F[s, i, j])

        return unit, _i, _j
