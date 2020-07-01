class File_IO():

    '''文件IO'''

    def __init__(self, root, handle):
        self.__M = handle
        self.root = root

    def input_gaussian_file(self, path):
        inp = open(path, 'r')

        i = 0
        inp_f = {0:[]}
        for line in inp:
            l = line.split()
            if l == []:
                i += 1
                inp_f[i] = []
                continue
            inp_f[i].append(l)
        
        inp.close()
        
        self.__M.set_charge(inp_f[2][0][0])
        self.__M.set_multiplicity(inp_f[2][0][0])
        atoms = self.__M.modify_atoms()
        for l in inp_f[2][1:]:
            atoms.append(l[0:1] + [float(ll) for ll in l[1:]])

        bonds = self.__M.modify_bonds()
        [bonds.append(dict()) for i in range(len(atoms))]
        for i in range(len(inp_f[3])):
            l = inp_f[3][i]
            for j in range((len(l)-1)//2):
                bonds[i][int(l[2*j+1])-1] = float(l[2*j+2])
                bonds[int(l[2*j+1])-1][i] = float(l[2*j+2])
    
    def output_gaussian_file(self, path):
        out = open(path, 'w')
        out.writelines('\n\n0 1\n')

        for a in self.__M.modify_atoms():
            out.write(' {:<15s}'.format(a[0]) + ' ' + ' '.join(['{:>13s}'.format('{:+.8f}'.format(s)) for s in a[1:]]) + '\n')
        
        out.writelines('\n')

        for i, b in enumerate(self.__M.modify_bonds()):
            out.write(' %d '%(i+1) + ' '.join([str(j+1) + ' ' + '{:.1f}'.format(b[j]) for j in b if i<j]) + '\n')
        
        out.writelines('\n')

        out.close()