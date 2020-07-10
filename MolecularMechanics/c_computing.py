import numpy as np
from ctypes import *

lib = None
fortran = 1
path = 'MolecularMechanics/'

def init(bit):
    global lib
    if bit == 32:
        lib = WinDLL(path + 'mm_c.dll')
    else:
        if fortran:
            lib = CDLL(path + 'mm_f_64.dll')
        else:
            lib = WinDLL(path + 'mm_c_64.dll')

def fdihedral_c(num2, array1, array2, array3, array4, dihedraldata):
    ans = np.asarray((0, 0, 0), dtype = 'double')
    p_ans = ans.ctypes.data_as(c_char_p)
    array1, array2, array3, array4, dd = map(lambda x: np.asarray(x).ctypes.data_as(c_char_p), (array1, array2, array3, array4, dihedraldata[1]))
    if fortran:
        lib.FDIHEDRAL_C(array1, array2, array3, array4, dd, byref(c_long(len(dihedraldata[1]))), p_ans)
    else:
        lib.fdihedral_c(array1, array2, array3, array4, dd, len(dihedraldata[1]), p_ans)
    return ans

def fdihedral_s(num1, array1, array2, array3, array4, dihedraldata):
    ans = np.asarray((0, 0, 0), dtype = 'double')
    p_ans = ans.ctypes.data_as(c_char_p)
    array1, array2, array3, array4, dd = map(lambda x: np.asarray(x).ctypes.data_as(c_char_p), (array1, array2, array3, array4, dihedraldata[1]))
    if fortran:
        lib.FDIHEDRAL_S(array1, array2, array3, array4, dd, byref(c_long(len(dihedraldata[1]))), p_ans)
    else:
        lib.fdihedral_s(array1, array2, array3, array4, dd, len(dihedraldata[1]), p_ans)
    return ans