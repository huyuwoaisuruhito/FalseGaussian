"""
原子种类判断
"""


class AtomTypeError(TypeError):
    pass

class Atoms:
    atomtypelist = []
    
    def __init__(self, num, bondingmap):
        self.type = None
        sites = [i[0] for i in bondingmap.sites]
        name = sites[num]
        bondings = bondingmap.bondings
        
        v_membered = bondingmap.v_membered_rings
        vi_membered = bondingmap.vi_membered_rings
        
        if name == "C":
            self.type = Atoms.__carbon(num, bondings, sites, v_membered, vi_membered)
        elif name == "H":
            self.type = Atoms.__hydrogen(num, bondings, sites)
        elif name == "N":
            self.type = Atoms.__nitrogen(num, bondings, sites, v_membered, vi_membered)
        elif name == "O":
            self.type = Atoms.__oxygen(num, bondings, sites)
        elif name == "S":
            self.type = Atoms.__sulfur(num, bondings, sites)
        elif name in ("P", "F", "Cl", "Br", "I"):
            self.type = name
        else:
            raise AtomTypeError("sprry! there is no bonding data for this atom: {0}".format(name))
        
        try:
            Atoms.atomtypelist[num] = self.type
        except:
            Atoms.atomtypelist += [None]*len(bondingmap.bondings)
            Atoms.atomtypelist[num] = self.type
    
    @classmethod
    def gettypelist(cls):
        return cls.atomtypelist
    
    @staticmethod
    def __carbon(num, bondings, sites, v_membered, vi_membered):
        if len(bondings[num]) == 4:
            return "CT"
        elif len(bondings[num]) == 2:
            raise AtomTypeError(
                "sorry! there is no bonding data for this sp hybridized carbon: {0}".format(str(num)))
        elif len(bondings[num]) == 3:
            bonds = bondings[num]
            neibers = "".join([sites[i[0]] for i in bonds])
            
            in_v = False
            in_vi = False
            for ring in v_membered:
                if num in ring:
                    in_v = True
                    i = ring.index(num)
                    ringneibers = sites[ring[(i+1)%5]] + sites[ring[(i+4)%5]]
                    ringneiber_nums = (ring[(i+1)%5], ring[(i+4)%5])
                    ringneiber_neibers = "".join([
                        sites[j[0]] for k in ringneiber_nums
                        for j in bondings[k] if sites[k] == "N" and j[0] != num
                    ])
            for ring in vi_membered:
                if num in ring:
                    in_vi = True
                    i = ring.index(num)
                    vi_ringneibers = sites[ring[(i+1)%6]] + sites[ring[(i+5)%6]]
            
            aro = False
            if 1.5 in [i[1] for i in bonds]:
                aro = True
            
            if in_v and aro:
                if in_vi:
                    if neibers in "NCCNC":
                        return "CN"
                    else:
                        return "CB"
                else:
                    if ringneibers == "NN":
                        if "H" in ringneiber_neibers:
                            return "CR"
                        return "CK"
                    elif ringneibers in "CNC":
                        if "H" in ringneiber_neibers:
                            if "H" not in ringneibers:
                                return "CC"
                            return "CW"
                        return "CV"
                    return "C*"
            elif in_vi and aro:
                if vi_ringneibers == "NN":
                    return "CQ"
                return "CA"
            elif "O" in neibers or "N" in neibers:
                if neibers == "NNN":
                    return "CA"
                return "C"
            return "CM"
        raise AtomTypeError(
                "sorry! there is no bonding data for this carbon: {0}".format(str(num)))
        
    def __nitrogen(num, bondings, sites, v_membered, vi_membered):
        bonds = bondings[num]
        neibers = "".join([sites[i[0]] for i in bonds])
        neiber_nums = [i[0] for i in bonds]
        
        def is_sp3(bonds):
            for i in bonds:
                if i[1] != 1:
                    return False
            return True
        
        def is_cation(bonds): return sum([i[1] for i in bonds]) == 4
        
        if is_sp3(bonds):
            return "N3"
        if len(bonds) == 1:
            raise AtomTypeError(
                "sorry! there is no bonding data for this nitrogen: {0}".format(str(num)))
        
        in_v = False
        in_vi = False
        for ring in v_membered:
            if num in ring:
                in_v = True
        for ring in vi_membered:
            if num in ring:
                in_vi = True

        aro = False
        if 1.5 in [i[1] for i in bonds]:
            aro = True
        
        if in_v and aro:
            if len(bonds) == 3:
                if "H" in neibers:
                    return "NA"
                return "N*"
            return "NB"
        if in_vi and aro:
            return "NC"
        
        neibertypes = [Atoms.atomtypelist[i] for i in neiber_nums]
        if "C" in neibertypes:
            return "N"
        return "N2"
    
    def __hydrogen(num, bondings, sites):
        bonds = bondings[num]
        neiber_num = bonds[0][0]
        neiber = sites[neiber_num]
        
        if neiber == "N":
            if sum([i[1] for i in bondings[neiber_num]]) == 4:
                return "HP"
            return "H"
        if neiber == "O":
            return "HO"
        if neiber == "S":
            return "HS"
        if neiber == "C":
            carbontype = Atoms.atomtypelist[neiber_num]
            if carbontype == "CT":
                neiber_nums = [i[0] for i in bondings[neiber_num] if i[0] != num]
                ewg = 0
                for n in neiber_nums:
                    if sites[n] in ("P", "F", "Cl", "Br", "I", "N", "O", "S"):
                        ewg += 1
                        continue
                    if sites[n] == "C":
                        if Atoms.gettypelist()[n] in ("CQ", "C"):
                            ewg += 1
                            continue
                if ewg == 1:
                    return "H1"
                if ewg == 2:
                    return "H2"
                if ewg == 3:
                    return "H3"
                return "HC"
            if carbontype in ("CQ", "CK", "CR"):
                return "H5"
            if carbontype in ("CV", "CC", "CW", "C"):
                return "H4"
            return "HA"
        raise AtomTypeError(
                "sorry! there is no bonding data for this hydrogen: {0}".format(str(num)))
        
    def __oxygen(num, bondings, sites):
        bonds = bondings[num]
        neibers = "".join([sites[i[0]] for i in bonds])
        
        if sum([i[1] for i in bonds]) < 2:
            return "O2"
        if len(neibers) == 2:
            if "H" in neibers:
                return "OH"
            return "OS"
        return "O"
    
    def __sulfur(num, bondings, sites):
        bonds = bondings[num]
        neibers = "".join([sites[i[0]] for i in bonds])
        if "H" in neibers:
            return "SH"
        return "S"

