import ctypes
from timeit import default_timer as timer

import numpy as np
from scipy import special

Lib_boysGaussJacobi = ctypes.CDLL('lib/.build/boysGaussJacobi.lib')
Lib_boysGaussJacobi.boysGaussJacobi.argtype = (ctypes.c_int, ctypes.c_double)
Lib_boysGaussJacobi.boysGaussJacobi.restype = ctypes.c_double


def dot_product(R1, R2):
    return sum([(R1[i]-R2[i])**2 for i in range(3)])


def N(e, ang):
    N = (4*e)**(sum(ang))
    for i in range(3):
        N /= special.factorial2(2*ang[i]-1, exact=True)
    N *= ((2*e)/np.pi)**(1.5)
    N = N**(0.5)
    return N


def nuclear_repulsion(atoms):
    Vnn = 0.0
    for a, A in enumerate(atoms):
        for b, B in enumerate(atoms):
            if b > a:
                num = A[0] * B[0]
                den = np.sqrt(dot_product(A[1:], B[1:]))
                Vnn += num / den
    return Vnn


def buildS(basis, S):
    K = len(basis)
    for i in range(K):
        for j in range(i+1):
            Bi, Bj = basis[i], basis[j]
            for bi in Bi.cgf:
                for bj in Bj.cgf:
                    S[i][j] += bi.k * bj.k * \
                        overlap_integral(bi.exp, bj.exp, bi.R,
                                         bj.R, bi.ang, bj.ang)
            S[j][i] = S[i][j]
    return S


def buildT(basis, T):
    K = len(basis)
    for i in range(K):
        for j in range(i+1):
            Bi, Bj = basis[i], basis[j]
            for bi in Bi.cgf:
                for bj in Bj.cgf:
                    T[i][j] += bi.k * bj.k * \
                        kinetic_energy_integral(
                            bi.exp, bj.exp, bi.R, bj.R, bi.ang, bj.ang)
            T[j][i] = T[i][j]
    return T


def buildV(basis, atoms, V):
    K = len(basis)
    for i in range(K):
        for j in range(i+1):
            Bi, Bj = basis[i], basis[j]
            for bi in Bi.cgf:
                for bj in Bj.cgf:
                    for Z, _x, _y, _z in atoms:
                        R3 = np.array([_x, _y, _z])
                        V[i][j] += bi.k * bj.k * nuclear_attraction_integral(
                            bi.exp, bj.exp, bi.R, bj.R, bi.ang, bj.ang, R3, Z)
            V[j][i] = V[i][j]
    return V


def buildG(basis, G):
    N = 0
    K = len(basis)
    t = timer()
    for i in range(K):
        for j in range(K):
            for k in range(K):
                for l in range(K):
                    N += 1
                    if N % 1000 == 1:
                        t, dt = timer(), timer() - t
                        print('Computed ' + str(N-1) + ' of ' + str(K) + '^4=' +
                              str(K**4) + ' two-electron integrals in {:.4f} s.'.format(dt))
                    if G[i][j][k][l] != 0:
                        continue
                    Bi, Bj, Bk, Bl = basis[i], basis[j], basis[k], basis[l]
                    for bi in Bi.cgf:
                        for bj in Bj.cgf:
                            for bk in Bk.cgf:
                                for bl in Bl.cgf:
                                    Eri = electron_repulsion_integral(
                                        bi.exp, bj.exp, bk.exp, bl.exp, bi.R, bj.R, bk.R, bl.R, bi.ang, bj.ang, bk.ang, bl.ang)
                                    G[i][j][k][l] += bi.k * \
                                        bj.k * bk.k * bl.k * Eri
                    G[j][i][k][l] = G[i][j][k][l]
                    G[i][j][l][k] = G[i][j][k][l]
                    G[j][i][l][k] = G[i][j][k][l]
                    G[k][l][i][j] = G[i][j][k][l]  # 要求空间轨道为实函数
                    G[l][k][i][j] = G[i][j][k][l]  # 要求空间轨道为实函数
                    G[k][l][j][i] = G[i][j][k][l]  # 要求空间轨道为实函数
                    G[l][k][j][i] = G[i][j][k][l]  # 要求空间轨道为实函数
    return G


def buildG_P(basis, G):
    N = 0
    t = timer()

    def callback(res):
        nonlocal G, N, t
        if res:
            i, j, k, l = res[0]
            G[i][j][k][l] = res[1]
            G[j][i][k][l] = G[i][j][k][l]
            G[i][j][l][k] = G[i][j][k][l]
            G[j][i][l][k] = G[i][j][k][l]
            G[k][l][i][j] = G[i][j][k][l]  # 要求空间轨道为实函数
            G[l][k][i][j] = G[i][j][k][l]  # 要求空间轨道为实函数
            G[k][l][j][i] = G[i][j][k][l]  # 要求空间轨道为实函数
            G[l][k][j][i] = G[i][j][k][l]  # 要求空间轨道为实函数

        N += 1
        if N % 1000 == 1:
            t, dt = timer(), timer() - t
            print('Computed ' + str(N-1) + ' of ' + str(K) + '^4='
                  + str(K**4) + ' two-electron integrals. This cycle takes {:.4f} s.'.format(dt))

    K = len(basis)
    from multiprocessing import Pool
    p = Pool(8)
    for i in range(K):
        for j in range(K):
            for k in range(K):
                for l in range(K):
                    p.apply_async(__buildG_p, args=(
                        i, j, k, l, basis, G), callback=callback)
    p.close()
    p.join()
    return G


def __buildG_p(i, j, k, l, basis, G):
    if G[i][j][k][l] != 0:
        return None
    Bi, Bj, Bk, Bl = basis[i], basis[j], basis[k], basis[l]
    g = 0
    for bi in Bi.cgf:
        for bj in Bj.cgf:
            for bk in Bk.cgf:
                for bl in Bl.cgf:
                    Eri = electron_repulsion_integral(
                        bi.exp, bj.exp, bk.exp, bl.exp, bi.R, bj.R, bk.R, bl.R, bi.ang, bj.ang, bk.ang, bl.ang)
                    Eri *= bi.k * bj.k * bk.k * bl.k
                    g += Eri
    return ((i, j, k, l), g)


def overlap_integral(e1, e2, R1, R2, ang1, ang2):
    """
    [1|2]
    """
    S = __overlap_integral(e1, e2, R1, R2)
    S *= __S(e1, e2, R1, R2, ang1, ang2)
    S *= N(e1, ang1) * N(e2, ang2)
    return S


def kinetic_energy_integral(e1, e2, R1, R2, ang1, ang2):
    """
    [1|/Delta^{2}|2]
    """
    K = e2*(2*sum(ang2)+3) * __S(e1, e2, R1, R2, ang1, ang2)
    ang2 = ang2.copy()
    for i in range(3):
        ang2[i] += 2
        K -= 2*e2**(2) * __S(e1, e2, R1, R2, ang1, ang2)
        ang2[i] -= 2
    for i in range(3):
        ang2[i] -= 2
        K -= 1/2 * ang2[i]*(ang2[i]-1) * __S(e1, e2, R1, R2, ang1, ang2)
        ang2[i] += 2
    K *= N(e1, ang1) * N(e2, ang2)
    return K * __overlap_integral(e1, e2, R1, R2)


def nuclear_attraction_integral(e1, e2, R1, R2, ang1, ang2, R3, Z):
    """
    [1|V(3)|2]
    """
    P = (e1*R1+e2*R2)/(e1+e2)

    MPR3 = dot_product(P, R3)

    V = 0.0
    for l in range(ang1[0]+ang2[0]+1):
        for r in range(int(l/2)+1):
            for i in range(int((l-2*r)/2)+1):
                vx = __vi(l, r, i, e1+e2, R1[0], R2[0],
                          ang1[0], ang2[0], R3[0], P[0])

                for m in range(ang1[1]+ang2[1]+1):
                    for s in range(int(m/2)+1):
                        for j in range(int((m-2*s)/2)+1):
                            vy = __vi(
                                m, s, j, e1+e2, R1[1], R2[1], ang1[1], ang2[1], R3[1], P[1])

                            for n in range(ang1[2]+ang2[2]+1):
                                for t in range(int(n/2)+1):
                                    for k in range(int((n-2*t)/2)+1):
                                        vz = __vi(
                                            n, t, k, e1+e2, R1[2], R2[2], ang1[2], ang2[2], R3[2], P[2])
                                        nu = l+m+n-2*(r+s+t)-(i+j+k)
                                        F = __BoysFunction(nu, MPR3*(e1+e2))
                                        V += vx * vy * vz * F
    V *= -Z * 2*np.pi/(e1+e2) * __overlap_integral(e1, e2, R1, R2)
    V *= N(e1, ang1) * N(e2, ang2)
    return V


def electron_repulsion_integral(e1, e2, e3, e4, R1, R2, R3, R4, ang1, ang2, ang3, ang4):
    """
    [12|34]
    """
    ep = e1 + e2
    eq = e3 + e4

    delta = 1/ep + 1/eq

    P = (e1*R1+e2*R2)/(e1+e2)
    Q = (e3*R3+e4*R4)/(e3+e4)

    MPQ = dot_product(P, Q)

    G = 0.0
    for l in range(0, ang1[0]+ang2[0]+1):
        for r in range(0, int(l/2)+1):
            for lp in range(0, ang3[0]+ang4[0]+1):
                for rp in range(0, int(lp/2)+1):
                    for i in range(0, int((l+lp-2*r-2*rp)/2)+1):
                        gx = __gi(l, lp, r, rp, i, ang1[0], ang2[0], R1[0], R2[0],
                                  P[0], ep, ang3[0], ang4[0], R3[0], R4[0], Q[0], eq)

                        for m in range(0, ang1[1]+ang2[1]+1):
                            for s in range(0, int(m/2)+1):
                                for mp in range(0, ang3[1]+ang4[1]+1):
                                    for sp in range(0, int(mp/2)+1):
                                        for j in range(0, int((m+mp-2*s-2*sp)/2)+1):
                                            gy = __gi(
                                                m, mp, s, sp, j, ang1[1], ang2[1], R1[1], R2[1], P[1], ep, ang3[1], ang4[1], R3[1], R4[1], Q[1], eq)

                                            for n in range(0, ang1[2]+ang2[2]+1):
                                                for t in range(0, int(n/2)+1):
                                                    for np_ in range(0, ang3[2]+ang4[2]+1):
                                                        for tp in range(0, int(np_/2)+1):
                                                            for k in range(0, int((n+np_-2*t-2*tp)/2)+1):
                                                                gz = __gi(
                                                                    n, np_, t, tp, k, ang1[2], ang2[2], R1[2], R2[2], P[2], ep, ang3[2], ang4[2], R3[2], R4[2], Q[2], eq)

                                                                nu = l+lp+m+mp+n+np_-2 * \
                                                                    (r+rp+s+sp +
                                                                     t+tp)-(i+j+k)
                                                                F = __BoysFunction(
                                                                    nu, MPQ/delta)
                                                                G += gx * gy * gz * F
    G *= (2*np.pi**2)/(ep*eq)
    G *= np.sqrt(np.pi/(ep+eq))
    G *= __overlap_integral(e1, e2, R1, R2)
    G *= __overlap_integral(e3, e4, R3, R4)
    G *= N(e1, ang1) * N(e2, ang2) * N(e3, ang3) * N(e4, ang4)
    return G


def __overlap_integral(e1, e2, R1, R2):
    I = np.exp(-e1*e2/(e1+e2) * dot_product(R1, R2))
    return I


def __S(e1, e2, R1, R2, ang1, ang2):
    S = (np.pi/(e1+e2))**(3/2)
    P = (e1*R1+e2*R2)/(e1+e2)
    for i in range(3):
        s = 0
        for k in range(0, int((ang1[i]+ang2[i])/2)+1):
            s += __ck(2*k, ang1[i], ang2[i], P[i]-R1[i], P[i]-R2[i]) * \
                special.factorial2(2*k-1, exact=True) / (2*(e1+e2))**k
        S *= s
    return S


def __vi(a, b, c, g, r1, r2, l1, l2, r3, p):
    eps = 1/(4*g)
    vi = (-1)**a
    vi *= __ck(a, l1, l2, p-r1, p-r2)
    vi *= (-1)**c * special.factorial(a, exact=True)
    vi *= (p-r3)**(a-2*b-2*c) * eps**(b+c)
    vi /= special.factorial(b, exact=True)
    vi /= special.factorial(c, exact=True)
    vi /= special.factorial(a-2*b-2*c, exact=True)
    return vi


def __gi(l, lp, r, rp, i, lA, lB, Ai, Bi, Pi, gP, lC, lD, Ci, Di, Qi, gQ):
    delta = 1/(4*gP) + 1/(4*gQ)
    gi = (-1)**l
    gi *= __theta(l, lA, lB, Pi-Ai, Pi-Bi, r, gP) * \
        __theta(lp, lC, lD, Qi-Ci, Qi-Di, rp, gQ)
    gi *= (-1)**i * (2*delta)**(2*(r+rp))
    gi *= special.factorial(l+lp-2*r-2*rp, exact=True) * delta**i
    gi *= (Pi-Qi)**(l+lp-2*(r+rp+i))
    gi /= (4*delta)**(l+lp) * special.factorial(i, exact=True)
    gi /= special.factorial(l+lp-2*(r+rp+i), exact=True)
    return gi


def __theta(l, l1, l2, a, b, r, g):
    the = __ck(l, l1, l2, a, b) * special.factorial(l, exact=True) * g**(r-l)
    the /= special.factorial(r, exact=True)
    the /= special.factorial(l-2*r, exact=True)
    return the


def __ck(k, l, m, a, b):
    c = 0.0
    for i in range(l+1):
        for j in range(m+1):
            if j + i == k:
                c += special.binom(l, i) * special.binom(m,
                                                         j) * a**(l-i) * b**(m-j)
    return c


def __BoysFunction(nu, x):
    if x < 1e-6:
        return (2*nu+1)**(-1) - x*(2*nu+3)**(-1)
    else:
        return Lib_boysGaussJacobi.boysGaussJacobi(ctypes.c_int(nu), ctypes.c_double(x))
