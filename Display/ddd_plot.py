import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

F_RADII = {
    'H': 0.1, 'B': 0.15, 'C': 0.2, 'N': 0.2, 'O': 0.2, 'F': 0.2, 'Si': 0.3, 'P': 0.3, 'S': 0.3, 'Cl': 0.3,
}
F_COLOR = {
    'H': 'whitesmoke', 'B': 'darkgray', 'C': 'dimgray', 'N': 'blue', 'O': 'red', 'F': 'greenyellow', 'Si': 'whitesmoke', 'P': 'magenta', 'S': 'gold', 'Cl': 'limegreen',
}


class DDD_plot():

    '''matplotlib实现3d显示'''

    def __init__(self):
        self.fig = plt.figure(facecolor = 'mediumpurple')
        self.__selectable = []
    
    def init(self, molecule):
        self.ax = Axes3D(self.fig, facecolor=(0.5, 0.5, 0.797))
        self.plot(molecule)
        plt.ion()
    
    def plot(self, molecule):
        self.ax.set_axis_off()
        atoms = molecule.get_atoms()
        bonds = molecule.get_bonds()

        self.plot_atoms(atoms)

        self.ax.quiver([0, 0, 0], [0, 0, 0], [0, 0, 0],  [1, 0, 0], [0, 1, 0], [0, 0, 1])
        self.ax.text(1, 0, 0, "x")
        self.ax.text(0, 1, 0, "y")
        self.ax.text(0, 0, 1, "z")

        _max =max([max([abs(ll) for ll in l]) for l in zip(self.ax.get_xlim3d(), self.ax.get_ylim3d(), self.ax.get_zlim3d())]+[3])
        self.ax.set_xlim3d(-_max, _max)
        self.ax.set_ylim3d(-_max, _max)
        self.ax.set_zlim3d(-_max, _max)

        self.plot_bonds(atoms, bonds, 5/_max) 
    
    def re_plot(self, molecule):
        self.ax.clear()
        self.plot(molecule)
        plt.draw()
    
    def change_view_pos(self, e):
        D = 5
        if e.keysym == 'Up':
            self.ax.view_init(self.ax.elev + D, self.ax.azim)
        elif e.keysym == 'Down':
            self.ax.view_init(self.ax.elev - D, self.ax.azim)
        elif e.keysym == 'Left':
            self.ax.view_init(self.ax.elev, self.ax.azim - D)
        elif e.keysym == 'Right':
            self.ax.view_init(self.ax.elev, self.ax.azim + D)
        plt.draw()
    
    def change_view_dist(self, e):
        self.ax.dist -= e.delta/240
        plt.draw()

    def heigh_light_atom(self, a):
        if isinstance(a, int):
            ind = a
            artist = self.__selectable[a]
        else:
            artist = a
            ind = [s[0] for s in self.__selectable].index(artist)
        if not self.__selectable[ind][2]:
            artist.set_facecolors('aqua')
            self.__selectable[ind][2] = True
            return ind
        else:
            self.clear_high_light()
            return -1

    def clear_high_light(self):
        for i, (artist, color, flag) in enumerate(self.__selectable):
            if flag:
                artist.set_facecolors(color)
                self.__selectable[i][2] =  False
        return -1

    def plot_atoms(self, atoms):
        __selectable = []
        for i, atom in enumerate(atoms):
            r = F_RADII.get(atom[0], 0.2)
            c = F_COLOR.get(atom[0], 'fuchsia')
            f = False

            u = np.linspace(0, 2 * np.pi, 10)
            v = np.linspace(0, np.pi, 10)
            _x = r * np.outer(np.cos(u), np.sin(v)) + atom[1]
            _y = r * np.outer(np.sin(u), np.sin(v)) + atom[2]
            _z = r * np.outer(np.ones(np.size(u)), np.cos(v)) + atom[3]

            if self.__selectable != [] and self.__selectable[i][2]:
                c = 'aqua'
                f = True
            __selectable.append([self.ax.plot_surface(_x, _y, _z, color=c, alpha=0.5, picker=1), c, f])
            self.ax.text(atom[1], atom[2], atom[3], '[%d]' %i, fontsize='small')
        self.__selectable = __selectable
    
    def plot_bonds(self, atoms, bonding, l):
        for i, bonds in enumerate(bonding):
            A = atoms[i]
            if bonds != {}:
                delta = [1, 0, 0]                
                for j, bl in bonds.items():
                    if j<i-1: continue
                    B = atoms[j]
                    c_B = F_COLOR.get(B[0], 'fuchsia')
                    c_A = F_COLOR.get(A[0], 'fuchsia')

                    if bl == 1:
                        theta = np.linspace(0, 1, 8)
                        x = (A[1]-B[1])*theta + B[1]
                        y = (A[2]-B[2])*theta + B[2]
                        z = (A[3]-B[3])*theta + B[3]

                        self.ax.plot(x[:5], y[:5], z[:5], color = c_B, linewidth = l)
                        self.ax.plot(x[4:], y[4:], z[4:], color = c_A, linewidth = l)
                    
                    elif bl == 2 or bl == 1.5:    
                        if len(bonds)<2 and len(bonding[j])<2:
                            delta = [0.05, 0, 0]
                        else:
                            if len(bonds)<2:
                                A, B = B, A
                                bonds = bonding[j]
                            n1 = np.array([A[1]-B[1], A[2]-B[2], A[3]-B[3]])
                            for j, _bl in bonds.items():
                                if atoms[j] != B:
                                    C = atoms[j]
                                    n2 = np.array([A[1]-C[1], A[2]-C[2], A[3]-C[3]])
                            fxl = np.cross(n2, n1)
                            if sum(fxl) == 0:
                                fxl = [0.5, 0.5, 0.5]
                            delta = np.cross(n1, fxl)
                            delta = [l/np.sqrt(sum(map(lambda a: a**2, delta))) for l in delta]
                            delta = [delta[i]*0.05 for i in range(3)]

                        theta = np.linspace(0, 1, 8)
                        x = (A[1]-B[1])*theta + B[1]
                        y = (A[2]-B[2])*theta + B[2]
                        z = (A[3]-B[3])*theta + B[3]

                        if bl == 2:
                            self.ax.plot(x[:5] + delta[0], y[:5] + delta[1], z[:5] + delta[2], 
                            color = c_B, linewidth = l)
                            self.ax.plot(x[4:] + delta[0], y[4:] + delta[1], z[4:] + delta[2], 
                            color = c_A, linewidth = l)

                            self.ax.plot(x[:5] - delta[0], y[:5] - delta[1], z[:5] - delta[2], 
                            color = c_B, linewidth = l)
                            self.ax.plot(x[4:] - delta[0], y[4:] - delta[1], z[4:] - delta[2], 
                            color = c_A, linewidth = l)

                        else:
                            self.ax.plot(x[:5] + delta[0], y[:5] + delta[1], z[:5] + delta[2], 
                            color = c_B, linewidth = l)
                            self.ax.plot(x[4:] + delta[0], y[4:] + delta[1], z[4:] + delta[2], 
                            color = c_A, linewidth = l)

                            for i in range(4):
                                if i<=2:
                                    c = c_B
                                else:
                                    c = c_A
                                self.ax.plot(x[2*i:2*i+2] - delta[0], y[2*i:2*i+2] - delta[1], z[2*i:2*i+2] - delta[2], color = c, linewidth = l)
                    
                    elif bl == 3:
                        theta = np.linspace(0, 1, 8)
                        x = (A[1]-B[1])*theta + B[1]
                        y = (A[2]-B[2])*theta + B[2]
                        z = (A[3]-B[3])*theta + B[3]

                        n0 = np.array([0.5, 0.5, 0.5])
                        n1 = np.array([A[1]-B[1], A[2]-B[2], A[3]-B[3]])
                        delta = np.cross(n1, n0)
                        su = np.sqrt(delta[0]**2 + delta[1]**2 + delta[2]**2)
                        delta = [delta[i]/su/7.5/np.sqrt(l) for i in range(3)]

                        self.ax.plot(x[:5], y[:5], z[:5], color = c_B, linewidth = l)
                        self.ax.plot(x[4:], y[4:], z[4:], color = c_A, linewidth = l)

                        self.ax.plot(x[:5] + delta[0], y[:5] + delta[1], z[:5] + delta[2], 
                        color = c_B, linewidth = l)
                        self.ax.plot(x[4:] + delta[0], y[4:] + delta[1], z[4:] + delta[2], 
                        color = c_A, linewidth = l)

                        self.ax.plot(x[:5] - delta[0], y[:5] - delta[1], z[:5] - delta[2], 
                        color = c_B, linewidth = l)
                        self.ax.plot(x[4:] - delta[0], y[4:] - delta[1], z[4:] - delta[2], 
                        color = c_A, linewidth = l)