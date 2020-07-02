import copy
import numpy as np

class Molecule():

    '''分子结构与控制'''

    def __init__(self):
        self.__charge = 0
        self.__multiplicity = 1
        self.__atoms = []
        self.__bonds = []

        self.__N_AOB = {}

    def copy(self):
        return copy.deepcopy(self)

    def set_charge(self, c):
        self.__charge = c
    
    def set_multiplicity(self, m):
        self.__multiplicity = m
    
    def get_charge(self, c):
        return self.__charge
    
    def get_multiplicity(self, m):
        return self.__multiplicity

    def get_atoms(self):
        return copy.deepcopy(self.__atoms)

    def get_bonds(self):
        return copy.deepcopy(self.__bonds)

    def modify_atoms(self):
        return self.__atoms

    def modify_bonds(self):
        return self.__bonds

    def get_bond_length(self, a, b):
        if a == b:
            return 0
        return np.sqrt(sum(map(lambda A, B:(A-B)**2, self.__atoms[a][1:], self.__atoms[b][1:])))

    def get_bond_angle(self, a, o, b, f=1):
        r_NN = np.array([0, 0, -1])
        if len(set([a, o, b]))!=3:
            return 0
        OA = np.array(self.__atoms[a][1:]) - np.array(self.__atoms[o][1:])
        OB = np.array(self.__atoms[b][1:]) - np.array(self.__atoms[o][1:])
        ON = np.cross(OA, OB)
        ANGLE = np.arccos( OA.dot(OB)/(np.sqrt(OA.dot(OA)) * np.sqrt(OB.dot(OB))) )

        if f==1: # 0 to 2pai
            if np.isnan(ANGLE):
                return 0
            elif r_NN.dot(ON)>=0:
                return ANGLE
            elif r_NN.dot(ON)<0:
                return 2*np.pi - ANGLE
            else:
                return
        elif f==0: #0 to pai
            if np.isnan(ANGLE):
                return 0
            else:
                return ANGLE


    def get_dihedral_angle(self, a, b, c, d):
        r_NN = np.array([0, 0, -1])
        if len(set([a, b, c, d]))!=4:
            return 0
        AB = np.array(self.__atoms[b][1:]) - np.array(self.__atoms[a][1:])
        BC = np.array(self.__atoms[c][1:]) - np.array(self.__atoms[b][1:])
        CD = np.array(self.__atoms[d][1:]) - np.array(self.__atoms[c][1:])
        N_ABC = np.cross(AB, BC)
        N_BCD = np.cross(CD, BC)
        NN = np.cross(N_BCD, N_ABC)
        ANGLE = np.arccos( N_ABC.dot(N_BCD)/(np.sqrt(N_ABC.dot(N_ABC)) * np.sqrt(N_BCD.dot(N_BCD))) )
        if np.isnan(ANGLE):
            return 0
        elif NN.dot(BC)>=0:
            return ANGLE
        elif NN.dot(BC)<0:
            return -ANGLE
        else:
            return

    def get_bond_level(self, a, b):
        return self.__bonds[a].get(b, 0)

    def modify_bond_level(self, a, b, bl):
        if a == b:
            return
        if bl != 0:
            self.__bonds[a][b] = bl
            self.__bonds[b][a] = bl
        elif bl == 0:
            try:
                del self.__bonds[a][b]
                del self.__bonds[b][a]
            except KeyError:
                pass

    def __get_unconnected_atom(self, a, b):
        unc_atom = set()

        def __(a, b, root, past):
            for i in self.__bonds[b].keys():
                if i in past: continue
                if i == a:
                    root = set()
                    return 0
                else:
                    if not __(a, i, root, past + [i]): return 0
                    root.add(i)
            return 1

        for i in self.__bonds[b].keys():
            if i == a: continue
            root = set([i])
            if __(a, i, root, [b, i]):
                unc_atom = set.union(unc_atom, root)
        
        return list(unc_atom) + [b]

    def __get_ritation_matrix(self, v0, v1, theta):
        a, b, c = v0
        u, v, w = v1
        l = np.sqrt(sum(map(lambda x: x**2, (u, v, w))))
        u, v, w = map(lambda x: x/l, (u, v, w))

        uu, uv, uw = u * u, u * v, u * w
        vv, vw, ww = v * v, v * w, w * w
        au, av, aw = a * u, a * v, a * w
        bu, bv, bw = b * u, b * v, b * w
        cu, cv, cw = c * u, c * v, c * w

        costheta = np.cos(theta)
        sintheta = np.sin(theta)

        m = [[0 for i in range(4)]for i in range(4)]

        m[0][0] = uu + (vv + ww) * costheta
        m[0][1] = uv * (1 - costheta) + w * sintheta
        m[0][2] = uw * (1 - costheta) - v * sintheta
        m[0][3] = 0

        m[1][0] = uv * (1 - costheta) - w * sintheta
        m[1][1] = vv + (uu + ww) * costheta
        m[1][2] = vw * (1 - costheta) + u * sintheta
        m[1][3] = 0

        m[2][0] = uw * (1 - costheta) + v * sintheta
        m[2][1] = vw * (1 - costheta) - u * sintheta
        m[2][2] = ww + (uu + vv) * costheta
        m[2][3] = 0

        m[3][0] = (a * (vv + ww) - u * (bv + cw)) * (1 - costheta) + (bw - cv) * sintheta
        m[3][1] = (b * (uu + ww) - v * (au + cw)) * (1 - costheta) + (cu - aw) * sintheta
        m[3][2] = (c * (uu + vv) - w * (au + bv)) * (1 - costheta) + (av - bu) * sintheta
        m[3][3] = 1
        
        return np.matrix(m)

    def __modify_bond_length(self, a, b, f, **kw):
        if a == b:
            return
        l = kw['x']

        if f:
            b_group = self.__get_unconnected_atom(a, b)
        else:
            b_group = (b,)

        for b_g in b_group:
            BL = self.get_bond_length(a, b)
            d_AB = list(map(lambda A, B: (B-A)*(l-BL)/BL, self.__atoms[a][1:], self.__atoms[b][1:]))
            self.__atoms[b_g] = self.__atoms[b_g][0:1] + list(map(sum, zip(self.__atoms[b_g][1:], d_AB)))
    
    def modify_bond_length_A(self, a, b, **kw):
        self.__modify_bond_length(a, b, 0, **kw)

    def modify_bond_length_G(self, a, b, **kw):
        self.__modify_bond_length(a, b, 1, **kw)

    def __modify_bond_angle(self, a, o, b, f, **kw):
        if len(set([a, o, b]))!=3:
            return
        delta_angle = kw['delta'] * kw['mul']
        OA = np.array(self.__atoms[a][1:]) - np.array(self.__atoms[o][1:])
        OB = np.array(self.__atoms[b][1:]) - np.array(self.__atoms[o][1:])
        ON = np.cross(OA, OB)
        M = self.__get_ritation_matrix(self.__atoms[o][1:], ON, delta_angle)

        if f:
            b_group = self.__get_unconnected_atom(a, b)
        else:
            b_group = (b,)
        
        for b_g in b_group:
            pos = np.matrix(self.__atoms[b_g][1:]+[1]) * M
            self.__atoms[b_g] = self.__atoms[b_g][0:1] + list(pos.A[0][:-1])
    
    def modify_bond_angle_A(self, a, o, b, **kw):
        self.__modify_bond_angle(a, o, b, 0, **kw)
    
    def modify_bond_angle_G(self, a, o, b, **kw):
        self.__modify_bond_angle(a, o, b, 1, **kw)

    def __modify_dihedral_angle(self, a, b, c, d, f, **kw):
        if len(set([a, b, c, d]))!=4:
            return
        delta_angle = kw['delta']

        while delta_angle>np.pi or delta_angle<-np.pi:
            if delta_angle>0:
                delta_angle -= 2*np.pi
            elif delta_angle<0:
                delta_angle += 2*np.pi
        
        delta_angle *= kw['dmul']
        BC = np.array(self.__atoms[c][1:]) - np.array(self.__atoms[b][1:])
        M = self.__get_ritation_matrix(self.__atoms[b][1:], BC, delta_angle)

        if f:
            if kw['dmul']>0:
                b_group = self.__get_unconnected_atom(c, b)
            else:
                b_group = self.__get_unconnected_atom(b, c)
        else:
            b_group = (d,)

        for b_g in b_group:
            pos = np.matrix(self.__atoms[b_g][1:]+[1]) * M
            self.__atoms[b_g] = self.__atoms[b_g][0:1] + list(pos.A[0][:-1])
    
    def modify_dihedral_angle_A(self, a, b, c, d, **kw):
        self.__modify_dihedral_angle(a, b, c, d, 0, **kw)
    
    def modify_dihedral_angle_G(self, a, b, c, d, **kw):
        self.__modify_dihedral_angle(a, b, c, d, 1, **kw)


    def auto_set_O(self):
        delta = [0, 0, 0]
        sums = 0
        for n, x, y, z in self.__atoms:
            sums += n
            delta[0] += x * n
            delta[1] += y * n
            delta[2] += z * n
        if sums != 0:
            self.set_O(*[-d/sums for d in delta])
        else:
            return

    def set_O(self, x, y, z):
        self.__atoms = [ [atom[0], atom[1]+x, atom[2]+y, atom[3]+z] for atom in self.__atoms ]

    def add_atom(self, a):
        pass

    def add_bond(self, a, b):
        pass
