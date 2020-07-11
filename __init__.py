import Display.gui as gui
from HartreeFock.hfArchive import HFArchive
import HartreeFock.ci as ci

if __name__ == "__main__":
    # atoms = [[8, 0.000000, 0.000000, 0.227000], [1, 0.000000, 1.353000,-0.908000], [1, 0.000000,-1.353000,-0.908000]]
    atoms = [[1, 0.000000, 0.000000, 0.370065], [1, 0.000000, 0.000000,-0.370065]]
    hfa = HFArchive(2, atoms, '6-31g', 'H2', debug=0)
    hfa.init_molecular_integrals()
    hfa.RHF()
    # hfa.dump()
    # hfa = HFArchive.load('.temp/H2/H2.6-31g')
    ci.CI_A(hfa, 0)
    # gui.Main_windows()