import numpy as np
from scipy.linalg import schur

def aux_schur(H):
    # due to column major in armadillo memptr
    # print(np.transpose(H))
    [T,Z,ndim] = schur(np.transpose(H),output='complex',sort='rhp')
    Tr = np.real(T)
    Ti = np.imag(T)
    Zr = np.real(Z)    
    Zi = np.imag(Z)
    # print('T(python) : {}'.format(T))
    # print('Z(python) : {}'.format(Z))
    return [Tr, Ti, Zr, Zi]