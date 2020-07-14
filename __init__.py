import Display.gui as gui
import HartreeFock.hfArchive as hfa

HFArchive = hfa.HFArchive
A0 = 0.52917721067


def draw(basis, fresh):
    name = './.temp/Plot/plot.'+basis
    colors = ['#1f77b4',
              '#ff7f0e',
              '#2ca02c',
              '#d62728',
              '#9467bd',
              '#8c564b',
              '#e377c2',
              '#7f7f7f',
              '#bcbd22',
              '#17becf',
              '#1a55FF']

    import os.path
    import numpy as np
    
    L = np.linspace(0.35, 1.0, 30)
    L = np.append(L, np.linspace(1.01, 1.51, 50))
    L = np.append(L, np.linspace(1.52, 3.0, 20))
    Zero = np.zeros(L.shape)
    Erhf = np.zeros(L.shape)
    Euhf = np.zeros(L.shape)
    Efci = np.zeros(L.shape)
    # Eufci = np.zeros(L.shape)

    if (os.path.exists(name+'.R.npy') and
        os.path.exists(name+'.U.npy') and
        os.path.exists(name+'.UFCI.npy') and
            not fresh):
        Erhf = np.load(name+'.R.npy')
        Euhf = np.load(name+'.U.npy')
        Efci = np.load(name+'.FCI.npy')
        # Eufci = np.load(name+'.UFCI.npy')
    else:
        for i, l in enumerate(L):
            print(l)
            atoms = [[1, 0.000000, 0.000000, l/2],
                     [1, 0.000000, 0.000000, -l/2]]
            hfA = HFArchive(2, atoms, basis, 'H2', fresh=1)
            hfA.init_molecular_integrals()
            Erhf[i] = hfA.RHF()
            Efci[i] = hfA.CI(0)
            Euhf[i] = hfA.UHF()
            # Eufci[i] = hfA.CI(0)

        hfa.dump_matrix(name, 'R', Erhf)
        hfa.dump_matrix(name, 'U', Euhf)
        hfa.dump_matrix(name, 'FCI', Efci)
        # hfa.dump_matrix(name, 'UFCI', Eufci)

    # from scipy import interpolate
    # Erhf_smooth = interpolate.interp1d(L, Erhf, kind='cubic')
    # Euhf_smooth = interpolate.interp1d(L, Euhf, kind='cubic')
    # Efci_smooth = interpolate.interp1d(L, Efci, kind='cubic')
    # Eufci_smooth = interpolate.interp1d(L, Eufci, kind='cubic')
    # nL = np.linspace(L.min(), L.max(), 100)

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    fig.suptitle('$\mathrm{H}_2$ dissociation (6-31g)')

    ax.plot(L/A0, Zero, label="Zero", color=colors[0])

    # ax.plot(L/A0, Erhf+1, 'o', color=colors[1])
    # ax.plot(L/A0, Euhf+1, 'o', color=colors[2])
    # ax.plot(L/A0, Efci+1, 'o', color=colors[3])
    # ax.plot(L/A0, Eufci+1, 'o', color=colors[4])

    # ax.plot(nL/A0, Erhf_smooth(nL)+1, label="RHF", color=colors[1])
    # ax.plot(nL/A0, Euhf_smooth(nL)+1, label="UHF", color=colors[2])
    # ax.plot(nL/A0, Efci_smooth(nL)+1, label="FCI", color=colors[3])
    # ax.plot(nL/A0, Eufci_smooth(nL)+1, label="UFCI", color=colors[4])
    
    ax.plot(L/A0, Erhf+1, '-', label="RHF" ,color=colors[1])
    ax.plot(L/A0, Euhf+1, '-', label="UHF", color=colors[2])
    ax.plot(L/A0, Efci+1, '-', label="FCI", color=colors[3])
    # ax.plot(L/A0, Eufci+1, '-', label="UFCI", color=colors[4])
    ax.legend()
    plt.show()


if __name__ == "__main__":
    # atoms1 = [[8, 0.000000, 0.000000, 0.227000], [1, 0.000000, 1.353000,-0.908000], [1, 0.000000,-1.353000,-0.908000]]
    # hfA = HFArchive(10, atoms1, '3-21g', 'H2O', fresh=0)
    atoms2 = [[1, 0.000000, 0.000000, 0.370065], [1, 0.000000, 0.000000,-0.370065]]
    # hfA = HFArchive(2, atoms2, '6-31g(d,p)', 'H2', fresh=1)
    hfA = HFArchive(2, atoms2, '3-21g', 'H2', fresh=1)
    # hfA = HFArchive(2, atoms2, 'sto-3g', 'H2', fresh=0)
    # atoms3 = [[1, 0.000000, 0.000000, 1], [1, 0.000000, 0.000000, -1]]
    # hfA = HFArchive(2, atoms3, 'sto-3g', 'H2', fresh=1)

    hfA.init_molecular_integrals()
    hfA.RHF()
    hfA.CI(0)
    # hfA.UHF()
    # hfA.CI(0)
    # hfA.dump()
    # gui.Main_windows()
    
    # draw('6-31g(d,p)', 0)