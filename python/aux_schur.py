import numpy as np
from scipy.linalg import schur

def aux_schur(_H):
    # due to column major in armadillo memptr
    # print(np.transpose(H))
    H = np.transpose(_H)
    T,Z,ndim = schur(H,output='complex',sort='lhp')
    Tr = np.real(T)
    Ti = np.imag(T)
    Zr = np.real(Z)    
    Zi = np.imag(Z)
    # print('T(python) : {}'.format(T))
    # print('Z(python) : {}'.format(Z))
    # dim = np.shape(H)[0]
    # half_dim = int(dim/2)
    # print('DIM : {}, {}'.format(half_dim, dim))
    # ISS = Z[:,0:half_dim]
    # X1 = ISS[0:half_dim,:]
    # X2 = ISS[half_dim:dim,:]
    # print('ISS : {}'.format(ISS))
    # print('X1 : {}'.format(X1))
    # print('X2 : {}'.format(X2))
    # print('H size : {}'.format(np.shape(H)))

    # X1_inv = np.linalg.inv(X1)
    # X = np.dot(X2, X1_inv)
    # temp1 = np.concatenate((X, -1*np.eye(half_dim)),axis=1)
    # temp2 = np.concatenate((np.eye(half_dim),X),axis=0)

    # print('Temp1 : {}'.format(temp1))
    # print('Temp2 : {}'.format(temp2))

    # H_c = H.astype(complex)
    # check1 = np.dot(H_c,temp2)
    # check = np.dot(temp1,check1)
    # print('CHECK : {}'.format(check))

    # H = np.matrix([[-3,2,0,0],
    #      [-2,1,0,-1],
    #      [0,0,3,2],
    #      [0,0,-2,-1]])
    # t,z,ndim = schur(H,output='complex',sort='lhp')
    # dim = np.shape(H)[0]
    # half_dim = int(dim/2)
    # iss = z[:,0:half_dim]
    # x1 = iss[0:half_dim,:]
    # x2 = iss[half_dim:dim,:]
    # x = np.dot(x2,np.linalg.inv(x1))
    # H_c = H.astype(complex)
    # temp1 = np.concatenate((x, -1*np.eye(half_dim)),axis=1)
    # temp2 = np.concatenate((np.eye(half_dim),x),axis=0)
    # check = np.dot(temp1,np.dot(H_c,temp2))
    # print('Unitary ? : {}'.format(np.dot(z,np.transpose(z))))
    # print('T(python) : {}'.format(t))
    # print('Z(python) : {}'.format(z))
    # print('ISS : {}'.format(iss))
    # print('X1 : {}'.format(x1))
    # print('X2 : {}'.format(x2))
    # print('H size : {}'.format(np.shape(H)))
    # print('CHECK : {}'.format(check))

    return [Tr, Ti, Zr, Zi]
