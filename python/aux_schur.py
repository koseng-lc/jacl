import numpy as np
from scipy.linalg import schur

def aux_schur(_H,_sort):
    # due to column major in armadillo memptr
    H = np.transpose(_H)
    T,Z,ndim = schur(H,output='complex',sort=_sort)
    Tr = np.real(T)
    Ti = np.imag(T)
    Zr = np.real(Z)    
    Zi = np.imag(Z)
    return [Tr, Ti, Zr, Zi]
