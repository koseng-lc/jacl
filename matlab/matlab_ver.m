% Nominal Parameters

% Viscous Friction Constant (N m s / rad)
Bm = 0.000048

% Rotor Inertia (kg m^2)
Jm = 0.00000072

% Torque Constant (N m / A)
Ki = 0.052

% Armature Inductance (H)
La = 0.00062

% Back EMF Constant(V s / rad)
Kb = 0.052

% Armature Resistance (Ohm)
Ra = 2.07

A = [0      1      0;
     0 -Bm/Jm  Ki/Jm;
     0 -Kb/La -Ra/La]

B = [   0     0;
        0 -1/Jm;
     1/La     0]

C = [1 0 0;
     0 1 0]
    
D = [0 0; 0 0]

t = 0:0.0001:0.01;
u1 = ones(size(t)) * 0.1
u2 = ones(size(t)) * -5
u = [u1' u2']

sys = ss(A,B,C,D)
% step(sys)
lsim(sys,u,t)