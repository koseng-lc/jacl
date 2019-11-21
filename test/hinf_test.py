import numpy as np
from scipy.signal import ss2tf

#-- Nominal Parameters

# Viscous Friction Constant (N m s / rad)
Bm = 0.000048

# Rotor Inertia (kg m^2)
Jm = 0.0000072

# Torque Constant (N m / A)
Ki = 0.052

# Armature Inductance (H)
La = 0.00062

# Back EMF Constant(V s / rad)
Kb = 0.052

# Armature Resistance (Ohm)
Ra = 2.07

#--

'''
    Input : |Voltage|
            |Torque |
'''

# -- LFT Matrix - P-Delta

A = np.matrix([
    [0,      1,      0],
    [0, -Bm/Jm,  Ki/Jm],
    [0, -Kb/La, -Ra/La]
])

B1 = np.matrix([
    [    0,     0,    0,     0,     0,     0],
    [-1/Jm, -1/Jm, 1/Jm,     0,     0,     0],
    [    0,     0,    0, -1/La, -1/La, -1/La]
])

B2 = np.matrix([
    [    0,     0],
    [    0, -1/Jm],
    [-1/La,     0]
])

C1 = np.matrix([
    [0,      1,      0],
    [0, -Bm/Jm,  Ki/Jm],
    [0,      0,      1],
    [0, -Kb/La, -Ra/La],
    [0,      1,      0],
    [0,      0,      1]
])

C2 = np.matrix([
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]
])

D11 = np.matrix([
    [    0,     0,    0,     0,     0,     0],
    [-1/Jm, -1/Jm, 1/Jm,     0,     0,     0],
    [    0,     0,    0,     0,     0,     0],
    [    0,     0,    0, -1/La, -1/La, -1/La],
    [    0,     0,    0,     0,     0,     0],
    [    0,     0,    0,     0,     0,     0],
])

D12 = np.matrix([
    [   0,     0],
    [   0, -1/Jm],
    [   0,     0],
    [1/La,     0],
    [   0,     0],
    [   0,     0]
])

D21 = np.zeros((3, 6))

D22 = np.zeros((3, 2))

#--

B = np.concatenate((B1, B2), axis = 1)

C = np.concatenate((C1, C2), axis = 0)

D = np.concatenate(
    (np.concatenate((D11, D12), axis = 1), np.concatenate((D21, D22), axis = 1)),
    axis = 0
)

def tf(_s, _A, _B, _C, _D):
    temp = (_s * np.eye(np.size(_A, 0)) - _A)
    return _C * np.linalg.inv(temp) * _B + _D

if __name__ == "__main__":
    test = tf(10, A, B, C, D)

    dc_motor_tf = ss2tf(A,B,C,D)

    print('A : {}'.format(A))
    print('B : {}'.format(B))
    print('C : {}'.format(C))
    print('D : {}'.format(D))

    print('TF(10) : {} '.format(test))