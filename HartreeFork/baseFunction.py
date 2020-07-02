import numpy as np

def doubleFact(x):
    ans=1
    for i in range(1,x+1):
        if i%2 == x%2:
            ans*=i
    return ans


class PGF:
    def __init__(self, k, exp, R, ang=None, ang_R=None):
        self.k = k
        self.exp = exp
        self.R = R

        if ang_R == None:
            self.ang_R = []
            for i, a in enumerate(ang):
                self.ang_R.append([] + [R[i]]*a)
            self._normalization()
        else:
            self.ang_R = ang_R


    def __mul__(self, b):
        a = self
        if isinstance(b, PGF):
            ang_R = [a_ang + b_ang for a_ang, b_ang in zip(a.ang_R, b.ang_R)]
            exp = a.exp + b.exp
            k = np.exp( -(a.exp*b.exp)/(a.exp+b.exp) * np.sum((a.R-b.R)**2) ) * a.k * b.k
            R = (a.exp*a.R + b.exp*b.R)/(a.exp+b.exp)
            return PGF(k, exp, R, ang_R=ang_R)
        elif isinstance(b, float):
            return PGF(b*a.k, a.exp, a.R, ang_R=a.ang_R)
        elif isinstance(b, int):
            return PGF(b*a.k, a.exp, a.R, ang_R=a.ang_R)
        else:
            return -1

    def _normalization(self):
        self.k /= np.sqrt((self * self).PGF_integration(0))

    def PGF_integration(self, f=1):
        ans = 1
        for r, a in zip(self.R, self.ang_R):
            if len(a)==0:
                ans *= np.sqrt(np.pi/self.exp)                                                                  #s
            if len(a)==1:
                ans *= np.sqrt(np.pi/self.exp) * (r-a[0])                                                       #p
            if len(a)==2:
                ans *= np.sqrt(np.pi/self.exp) * ( (r-a[0])*(r-a[1]) + 1/(2*self.exp) )                         #d
            if len(a)==3:
                ans *= np.sqrt(np.pi/self.exp) * ( (r-a[0])*(r-a[1])*(r-a[2]) + (3*r-sum(a))/(2*self.exp))      #f
        if f:
            return ans * self.k
        else:
            return ans


class CGF:
    def __init__(self, *arg, R=None, ang=None, ang_R=None, cgf=None, k=1):
        if cgf != None and ang_R != None:
            self.__init_2(k, ang_R, cgf)
        elif arg != None and tuple(R) != None and ang != None:
            self.__init_1(k, R, *arg, ang=ang)
        else:
            raise RuntimeError('CGF:Bad Argument')

    def __init_1(self, k, R, *arg, ang=None):
        self.ang_R = []
        for i, a in enumerate(ang):
            self.ang_R.append([] + [R[i]]*a)
        
        self.cgf = []
        for a in arg:
            self.cgf.append(PGF(*a, R, ang_R=self.ang_R))
        
        self.k = 1
        self._normalization()

    def __init_2(self, k, ang_R, cgf):
        self.ang_R = ang_R
        self.cgf = cgf
        self.k = k

    def __mul__(self, b):
        a = self
        if isinstance(b, CGF):
            ang_R = [a_ang + b_ang for a_ang, b_ang in zip(a.ang_R, b.ang_R)]
            n_cgf = []
            for a_pgf in a.cgf:
                for b_pgf in b.cgf:
                    n_cgf.append(a_pgf * b_pgf)
            return CGF(ang_R=ang_R, cgf=n_cgf, k=a.k*b.k)
        elif isinstance(b, float):
            return CGF(ang_R=ang_R, cgf=n_cgf, k=a.k*b)
        elif isinstance(b, int):
            return CGF(ang_R=ang_R, cgf=n_cgf, k=a.k*b)
        else:
            return -1

    def _normalization(self):
        self.k /= np.sqrt((self * self).CGF_integration(0))

    def CGF_integration(self, f=1):
        ans = 0
        for pgf in self.cgf:
            ans += pgf.PGF_integration()
        if f:
            return ans * self.k
        else:
            return ans


if __name__ == "__main__":
    p1 = PGF(1, 0.7, np.array((0,-1/2,0)), ang=np.array((0,0,0)))
    p2 = PGF(1, 1/3, np.array((0,1,0)), ang=np.array((0,0,0)))
    print((p1 * p1).PGF_integration())
    print((p2 * p2).PGF_integration())
    pp = p1 * p2
    print((pp * pp).PGF_integration())
