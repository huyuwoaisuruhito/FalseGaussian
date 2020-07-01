"""
用于计算的分子信息传输（原子种类及坐标，单向成键图，双向成键图，环系数据）
"""

class Bondingmaps:
    def __init__(self, molecule):            #"待改为molecules"
        self.sites = molecule.get_atoms()
        self.bondings = []
        self.relations = []
        self.bondinglist = []
        
        for a, bondsdict in enumerate(molecule.get_bonds()):
            self.bondings.append([])
            self.relations.append([])
            self.bondinglist.append([])
            for b, bondlevel in bondsdict.items():
                self.bondings[a].append([b, bondlevel])
                self.relations[a].append(b)
                if b > a:
                    self.bondinglist[a].append([b, bondlevel])
        
        self.rings = Bondingmaps.__findrings(self.relations)
        self.v_membered_rings = filter(lambda a: len(a)==5 ,self.rings)
        self.vi_membered_rings = filter(lambda a: len(a)==6 ,self.rings)
        
    def __getframework(relation):
        ends = []
        for i in range(len(relation)):
            if len(relation[i]) <= 1:
                ends.append(i)
        
        framework = []
        for i, bonds in enumerate(relation):
            if i in ends:
                framework.append([])
            else:
                framebonds = []
                for bond in bonds:
                    if bond not in ends:
                        framebonds.append(bond)
                framework.append(framebonds)
        
        return framework
    
    def __findrings(relations):
        ends = []
        for i in range(len(relations)):
            if len(relations[i]) <= 1:
                ends.append(i)
        
        framework = relations
        while True:
            another = Bondingmaps.__getframework(framework)
            if another == framework:
                break
            framework = another

        atoms = [0] * len(relations)
        crossings = [i for i, bonds in enumerate(framework) if len(bonds) >= 3]
        
        if len(crossings) == 0:
            if len(framework) != 0:
                rings = [[i for i, bonds in enumerate(framework) if len(bonds) != 0]]
                return rings
            return []
        
        edges = dict()
        for start in crossings:
            atoms[start] = 1
            for p in framework[start]:
                if atoms[p]:
                    continue
                
                path = []
                count = 0
                while p not in crossings:#death loop with nul
                    path.append(p)
                    atoms[p] =  1
                    count += 1
                    for at in framework[p]:
                        if at in crossings and count != 1:
                            p = at
                        if not atoms[at]:
                            p = at
                
                if start <= p:
                    key = (start, p)
                    if key not in edges:
                        edges[key] = []
                    edges[key].append(path)
        
        rings = []
        shortedges = dict()
        for key, paths in edges.items():
            if key[0] == key[1]:
                for path in paths:
                    rings.append([key[0]] + path)
            else:
                newpaths = sorted(paths, key=lambda path: len(path))
                while len(newpaths) != 1:
                    rings.append([key[0]] + newpaths[0] + [key[1]] + list(reversed(newpaths.pop())))
                shortedges[key] = newpaths[0]
        
        ways = []
        for i in range(len(relations)):
            ways.append([])
        for key in shortedges:
            a, b = key[0], key[1]
            ways[a].append(b)
            ways[b].append(a)
        
        noderings = Bondingmaps.__findrings(ways)
        for nodering in noderings:
            ring = []
            for i, node in enumerate(nodering):
                key = (nodering[i-1], node)
                if key[0] > key[1]:
                    key = (node, nodering[i-1])
                ring += shortedges[key]
                ring.append(node)
            if len(ring) != 0:
                rings.append(ring)
        
        return rings

