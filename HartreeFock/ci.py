from timeit import default_timer as timer

import numpy as np
from scipy import linalg


def CI_A(HFArchive, level):
    return CI(HFArchive.N, HFArchive.K,
              HFArchive.Hc, HFArchive.G,
              HFArchive.C, HFArchive.Vnn, 
              HFArchive.F, level,
              HFArchive.debug, HFArchive)


def CI(N, K, Hc, G, C, Vnn, F_, level, debug=0, HFArchive=None):
    '''
    CI方法，已测试
    level: 计算到第level重激发，为0时计算FCI
    '''
    t = timer()
    print('\n========== Begin CI ==========\n')

    if 2 == 1 and len(N) == 1:
        N = [N[0]//2, N[0]//2]
    elif 2 == 2:
        if len(N) == 1:
            N = [N[0]-N[0]//2, N[0]//2]
        elif len(N) == 2:
            pass
        else:
            raise RuntimeError()
    else:
        raise RuntimeError()

    exc_ob = [[], [], []]
    for s in range(2):

        if level == 0 or level > min(N[s], K-N[s]):
            max_c = min(N[s], K-N[s])
        else:
            max_c = level

        exc_ob_s = [[[], []]]
        __make_ci_c(s, N[s], K, max_c, n=N[s], res=exc_ob_s)
        exc_ob[s] = exc_ob_s[1:]
    for n in exc_ob[0]:
        for m in exc_ob[1]:
            zd = [n[0]+m[0], n[1]+m[1]]
            exc_ob[2].append(zd)

    exc_ob = [[[], []]] + exc_ob[0] + exc_ob[1] + exc_ob[2]
    exc_ob.sort(key=lambda x: len(x[0]))

    # exc_ob = [a for a in exc_ob if a[0][0]==0]
    if level:
        exc_ob = [_ for _ in exc_ob if len(_[0]) <= level]

    CI_K = len(exc_ob)

    print(exc_ob)
    t, dt = timer(), timer() - t
    print('Make {} excited determinant wavefunctions in {:.4f} s'.format(CI_K, dt))

    E = __naive_method(sum(N), N, K, Hc, G, F_, C,
                        level, CI_K, exc_ob, debug)

    print('Solved CI matrix in {:.4f} s'.format(dt))
    print('\n====== CI in level:{} ======'.format(level))
    print('\nE = Eel + Vnn')
    print('  = {:.6f} + {:.6f}'.format(E, Vnn))
    print('  = {:.6f} Hartrees\n'.format(E+Vnn))

    return E + Vnn


def __make_ci_c(s, N, K, max_c, count=0, m=0, n=0, a=[], r=[], res=[]):
    if count == max_c:
        return
    for i in range(m, N):
        for j in range(n, K):
            a.append(i*2+s)
            r.append(j*2+s)
            res.append([a.copy(), r.copy()])
            __make_ci_c(s, N, K, max_c, count+1, i+1, j+1, a, r, res)
            a.pop()
            r.pop()


def __naive_method(sN, N, K, Hc, G, F_, C, level, CI_K, exc_ob, debug, **kw):

    def __callback(res):
        nonlocal CI_M, n, t
        i, j = res[1:3]
        CI_M[i, j] = res[0]
        CI_M[j, i] = CI_M[i, j]
        n += 1
        if n % 1000 == 1:
            dt = timer() - t
            print('Computed ' + str(n-1) + ' of ' + str(CI_K) + '^2='
                  + str(CI_K**2) + ' CI matrix unit in {:.4f} s.'.format(dt))
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
                                        N, K, Hc, G, i, j, F_),
                                  callback=__callback)
        p.close()
        p.join()
    else:
        for i, exc_ob1 in enumerate(exc_ob):
            for j, exc_ob2 in enumerate(exc_ob):
                if j <= i:
                    CI_M[i, j] = __compute_ci_matrix_unit(
                        exc_ob1, exc_ob2, C, sN, N, K, Hc, G, i, j, F_)[0]
                    CI_M[j, i] = CI_M[i, j]

    t, dt = timer(), timer() - t
    print('Compute {}x{} CI matrix in {:.4f} s'.format(CI_K, CI_K, dt))
    print(CI_M)
    print()

    e, C = linalg.eigh(CI_M)

    print('\ne: \n', e)
    print('\nC:\n ', C)
    return e[0]


def __compute_ci_matrix_unit(exc_ob1, exc_ob2, C, sN, N, K, Hc, G, _i, _j, F_):
    '''
    本方法针对自旋轨道的CI，注意用UHF轨道计算的CI能量比用RHF轨道计算的CI能量要高
    '''
    if C.shape[0] == 1:
        C = np.append(C, np.copy(C), axis=0)
        F_ = np.append(F_, np.copy(F_), axis=0)

    pC1 = np.array([i for i in range(K*2)])
    for a, r in zip(*exc_ob1):
        pC1[[r, a]] = pC1[[a, r]]

    pC2 = np.array([i for i in range(K*2)])
    for a, r in zip(*exc_ob2):
        pC2[[r, a]] = pC2[[a, r]]

    dC = [(c1//2, c1 % 2, c2//2, c2 % 2)
          for c1, c2 in zip(pC1[:sN], pC2[:sN]) if c1 != c2]

    unit = 0
    if len(dC) > 2:
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
                        if sm == sn and sp == sq:
                            unit -= (C[sm, i, m] * C[sn, j, n] * C[sp, k, p]
                                     * C[sq, l, q] * G[i, j, k, l])  # [mn|pq]=K
        return unit, _i, _j

    P = np.zeros((2, K, K))
    for i in range(K):
        for j in range(K):
            for s in range(2):
                for a in range(N[s]):
                    na1 = pC1[2*a+s]//2
                    na2 = pC2[2*a+s]//2
                    if na1 == na2:
                        P[s, i, j] += C[s, i, na1] * C[s, j,na2]

    F = np.zeros((2, K, K))
    for s in range(2):
        for i in range(K):
            for j in range(K):
                F[s, i, j] = Hc[i, j]
                for k in range(K):
                    for l in range(K):
                            F[s, i, j] -= P[s, k, l] * G[i, k, j, l]
                            F[s, i, j] += P[0, k, l] * G[i, j, k, l]
                            F[s, i, j] += P[1, k, l] * G[i, j, k, l]

    if len(dC) == 1:
        m, sm, p, sp = dC[0]
        if sm != sp:
            return unit, _i, _j
        for i in range(K):
            for j in range(K):
                unit += C[sm, i, m] * C[sp, j, p] * F[sm, i, j]
        return unit, _i, _j

    if len(dC) == 0:
        for s in range(2):
            for i in range(K):
                for j in range(K):
                    unit += 1/2 * P[s, j, i] * (Hc[i, j] + F[s, i, j])
        return unit, _i, _j
