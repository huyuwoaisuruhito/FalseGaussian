import Display.gui as gui
from HartreeFock.hfArchive import HFArchive
import HartreeFock.ci2 as ci

if __name__ == "__main__":
    atoms1 = [[8, 0.000000, 0.000000, 0.227000], [1, 0.000000, 1.353000,-0.908000], [1, 0.000000,-1.353000,-0.908000]]
    hfa = HFArchive(10, atoms1, '3-21g', 'H2O', debug=0)
    # atoms2 = [[1, 0.000000, 0.000000, 0.370065], [1, 0.000000, 0.000000,-0.370065]]
    # hfa = HFArchive(2, atoms2, 'sto-3g', 'H2', debug=0)
    hfa.init_molecular_integrals()
    hfa.RHF()
    hfa.dump()
    # hfa = HFArchive.load('.temp/H2/H2.6-31g')
    hfa.CI(2)
    hfa.dump()
    # gui.Main_windows()