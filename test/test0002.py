'''
RHF & CI
'''

from HartreeFock.hfArchive import HFArchive

atoms = [[1, 0.000000, 0.000000, 0.370065], [1, 0.000000, 0.000000,-0.370065]]

hfA = HFArchive(2, atoms, '6-31g', 'H2', fresh=1)
hfA.init_molecular_integrals()

hfA.RHF()
hfA.CI(0)