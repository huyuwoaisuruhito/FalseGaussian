import timeit
from scipy import special

f_st = '''
import ctypes
import numpy as np
from scipy import special

Lib_boysGaussJacobi = ctypes.CDLL('lib/.build/boysGaussJacobi.lib')

def __ck(k, l, m, a, b):
    c = 0.0
    for i in range(l+1):
        for j in range(m+1):
            if j + i == k:
                c += special.binom(l, i) * special.binom(m, j) * a**(l-i) * b**(m-j)
    return c
def __theta(l, l1, l2, a, b, r, g):
    the = __ck(l, l1, l2, a, b) * special.factorial(l, exact=True) * g**(r-l)
    the /= special.factorial(r, exact=True)
    the /= special.factorial(l-2*r, exact=True)
    return the
'''

f_gi = '''
(l,lp,r,rp,i, lA,lB,Ai,Bi,Pi,gP, lC,lD,Ci,Di,Qi,gQ) = (0, 2, 0, 0, 1, 0, 0, 0.0, 0.0, 0.0, 370.4678, 1, 1, 0.0, 0.0, 0.0, 8.97914)
delta = 1/(4*gP) + 1/(4*gQ)
gi  = (-1)**l 
gi *= __theta(l,lA,lB,Pi-Ai,Pi-Bi,r,gP) * __theta(lp,lC,lD,Qi-Ci,Qi-Di,rp,gQ)
gi *= (-1)**i * (2*delta)**(2*(r+rp))
gi *= special.factorial(l+lp-2*r-2*rp,exact=True) * delta**i
gi *= (Pi-Qi)**(l+lp-2*(r+rp+i))
gi /= (4*delta)**(l+lp) * special.factorial(i,exact=True)
gi /= special.factorial(l+lp-2*(r+rp+i),exact=True)
'''

f_bf_1 = '''
def __BoysFunction(x, nu):
    if x < 1e-6:
        return (2*nu+1)**(-1) - x*(2*nu+3)**(-1)
    else:
        return Lib_boysGaussJacobi.boysGaussJacobi(ctypes.c_int(nu), ctypes.c_double(x))
__BoysFunction(2008.7598897159996, 1)
'''

f_bf_0 = '''
def __BoysFunction(x, nu):
    if x < 1e-6:
        return (2*nu+1)**(-1) - x*(2*nu+3)**(-1)
    else:
        return (1/2) * x**(-(nu+0.5)) * special.gamma(nu+0.5) * special.gammainc(nu+0.5,x)
__BoysFunction(2008.7598897159996, 1)
'''

f_bf_nf = '''
x = 2008.7598897159996
nu = 1
if x < 1e-6:
    (2*nu+1)**(-1) - x*(2*nu+3)**(-1)
else:
    (1/2) * x**(-(nu+0.5)) * special.gamma(nu+0.5) * special.gammainc(nu+0.5,x)
'''

factorial = '''
n=3
ans = 1
for i in range(n):
    ans *= i+1
'''


if __name__ == "__main__":
    # print(timeit.timeit(stmt=f_bf_0, setup=f_st, number=10000)) #0.0274 s
    # print(timeit.timeit(stmt=f_bf_1, setup=f_st, number=10000)) #0.0147 s
    # print(timeit.timeit(stmt=f_bf_nf, setup=f_st, number=10000)) #0.0289 s
    # print(timeit.timeit(stmt=f_gi, setup=f_st, number=10000))  # 0.5632 s
    print(timeit.timeit(stmt='special.factorial2(3, exact=True)', setup=f_st, number=90000))  # 0.1013 s
    print(timeit.timeit(stmt='special.factorial(3, exact=True)', setup=f_st, number=90000))  # 0.3542 s
    print(timeit.timeit(stmt=factorial, setup=f_st, number=90000))  # 0.0657 s
