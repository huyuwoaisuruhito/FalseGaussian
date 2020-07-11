import Display.gui as gui
import HartreeFock.RHF as rhf

if __name__ == "__main__":
    atoms = [[8, 0.000000, 0.000000, 0.227000], [1, 0.000000, 1.353000,-0.908000], [1, 0.000000,-1.353000,-0.908000]]
    # atoms = [[1, 0.000000, 0.000000, 0.370065], [1, 0.000000, 0.000000,-0.370065]]
    rhf.A_2_atom_unit(atoms)
    rhf.RHF(10, atoms, '3-21g', 'H2O')
    # gui.Main_windows()