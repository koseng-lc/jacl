'''
    author : koseng (Lintang)
    brief : Showing uncertainty impact on DC Motor
'''

import random
import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.linalg import expm

class DCMotorOpenLoop:
    def __init__(self, Bm, Jm, Ki, La, Kb, Ra):
        self.Bm = Bm
        self.Jm = Jm
        self.Ki = Ki
        self.La = La
        self.Kb = Kb
        self.Ra = Ra
        
        self.mod_Bm = self.Bm
        self.mod_Jm = self.Jm
        self.mod_Ki = self.Ki
        self.mod_La = self.La
        self.mod_Kb = self.Kb
        self.mod_Ra = self.Ra

    def perturbBm(self, delta):
        self.mod_Bm = self.Bm + delta

    def perturbJm(self, delta):
        self.mod_Jm = self.Jm + delta

    def perturbKi(self, delta):
        self.mod_Ki = self.Ki + delta
    
    def perturbLa(self, delta):
        self.mod_La = self.La + delta

    def perturbKb(self, delta):
        self.mod_Kb = self.Kb + delta

    def perturbRa(self, delta):
        self.mod_Ra = self.Ra + delta

    def unperturbAll(self):
        self.mod_Bm = self.Bm
        self.mod_Jm = self.Jm
        self.mod_Ki = self.Ki
        self.mod_La = self.La
        self.mod_Kb = self.Kb
        self.mod_Ra = self.Ra

    def A(self):
        return np.matrix([
            [0,                        1,                        0],
            [0, -self.mod_Bm/self.mod_Jm,  self.mod_Ki/self.mod_Jm],
            [0, -self.mod_Kb/self.mod_La, -self.mod_Ra/self.mod_La]
        ])

    def B(self):
        return np.matrix([
            [            0,              0],
            [            0, -1/self.mod_Jm],
            [1/self.mod_La,              0]
        ])

    def C(self):
        return np.matrix([
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]
        ])
    
    def D(self):
        return np.zeros((3,2))

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

def tf(_s, _A, _B, _C, _D):
    temp = (_s * np.eye(np.size(_A, 0)) - _A)
    return _C * np.linalg.inv(temp) * _B + _D

def bodePlot(title, _A, _B, _C, _D):
    freq = 2.0j * np.pi * np.logspace(0,5)
    # print('Freq : {}'.format(freq))
    magnitude = []
    phase = []
    for f in freq:
        val = tf(f, _A, _B, _C, _D)[0,0]
        phase_val = np.arctan2(np.imag(val), np.real(val)) * 180.0 / np.pi
        magnitude_val = 20 * np.log10(np.sqrt(np.imag(val) * np.imag(val) + np.real(val) * np.real(val)))        
        phase.append(phase_val)
        magnitude.append(magnitude_val)

        # Debug
        # print('Magnitude : {}'.format(magnitude_val))
        # print('Phase : {}'.format(phase_val))

    plt.clf()

    plt.subplot(211)
    plt.title(title)
    plt.xscale('log')
    # plt.yticks(np.arange(-150,150,step=25))
    plt.ylabel('Magnitude (dB)')
    plt.xlabel('Frequency (rad / sec)')
    plt.grid(True)
    plt.plot(np.imag(freq),magnitude)

    plt.subplot(212)
    plt.xscale('log')
    plt.ylabel('Phase (Degree)')
    plt.xlabel('Frequency (rad / sec)')
    plt.grid(True)
    plt.plot(np.imag(freq),phase)

    plt.show()

class LuenbergerObserver:
    def __init__(self):
        self.x_hat = 0
        self.y_hat = 0

def ssSimulation(title, _A, _B, _C, _D, _u):
    x = [0]
    y1 = [0]
    y2 = [0]

    plt.clf()
    plt.title(title)

    fig = plt.figure()
    
    sp1 = fig.add_subplot(211)
    sp2 = fig.add_subplot(212)

    line1, = sp1.plot(x,y1,color='b')
    line2, = sp2.plot(x,y2,color='r')

    fig.canvas.draw()
    fig.show()    

    largest_state_amp = 0.0
    largest_output_amp = 0.0

    sp1.set_ylim(bottom = -1.0, top = 1.0)
    sp2.set_ylim(bottom = -1.0, top = 1.0)

    # Time interval (seconds)
    d_time = 0.0001
    view_interval = 0.01

    max_data = view_interval / d_time

    # Time
    t = 0

    # State Transition Matrix
    state_trans = expm(_A * d_time)

    # Initial State
    init_state = np.zeros((_A.shape[0], 1))

    # State
    state = init_state
    last_state = state

    # Output
    output = np.zeros((_C.shape[1], 1))

    # Data
    state_data = [state[1]]
    output_data = [output[1]]
    time_data = [0]

    # Other stuff
    B_dt = _B * d_time
    total_state = 0
    total_output = 0

    # Input
    u = np.ones((_B.shape[1], 1))
    u[0] = 0.5 # Torque
    u[1] = -12 # Voltage
    last_u = np.zeros((_B.shape[1], 1))

    # Debug for validation
    print('A : {}'.format(_A))
    print('B : {}'.format(_B))
    print('C : {}'.format(_C))
    print('D : {}'.format(_D))

    print('State Transition : {}'.format(state_trans))
    print('Init State : {}'.format(init_state))
    print('Last State : {}'.format(last_state))
    print('State : {}'.format(state))
    print('B_dt : {}'.format(B_dt))
    print("Input : {}".format(u))
    print('Last Input : {}'.format(last_u))    

    while True:
        term1 = np.dot(state_trans, last_state)
        term2 = np.dot(B_dt, u)
        term3 = np.dot(state_trans, np.dot(B_dt, last_u))
        term4 = np.dot(state_trans, term2 + term3) * 0.5
        state = term1 + term4
        output = np.dot(_C, state) + np.dot(_D, u)
        # print('Term1 : {}'.format(term1))
        # print('Term2 : {}'.format(term2))
        # print('Term3 : {}'.format(term3))
        # print('Term4 : {}'.format(term4))
        print('State : {}'.format(state))
        print('Output : {}'.format(output))

        last_u = u
        last_state = state      

        total_state += np.abs(state[1])
        total_output += np.abs(output[1])

        state_data.append(state[1])
        output_data.append(output[1])
        time_data.append(t)

        line1.set_ydata(state_data)
        line1.set_xdata(time_data)

        line2.set_ydata(output_data)
        line2.set_xdata(time_data)  

        sp1.draw_artist(sp1.patch)
        sp1.draw_artist(line1)

        sp2.draw_artist(sp2.patch)
        sp2.draw_artist(line2)

        sp1.set_xlim(left = max(0, t - 0.75 * view_interval), right = max(view_interval, t + 0.25 * view_interval))
        sp2.set_xlim(left = max(0, t - 0.75 * view_interval), right = max(view_interval, t + 0.25 * view_interval))

        if(largest_output_amp < output[1]):
            largest_output_amp = output[1]

        if(largest_state_amp < state[1]):
            largest_state_amp = state[1]     

        output_avg = total_output / len(time_data)
        state_avg = total_state / len(time_data)

        sp1.set_ylim(bottom = -state_avg - state_avg * 0.75, top = state_avg + state_avg * 0.75)
        sp2.set_ylim(bottom = -output_avg - output_avg * 0.75, top = output_avg + output_avg * 0.75)

        fig.canvas.draw()
        fig.canvas.flush_events()

        t += d_time
        
        print('Time Data Size : {}'.format(len(time_data)))
        print('State Data Size : {}'.format(len(state_data)))
        print('Output Data Size : {}'.format(len(output_data)))

        print('Time : {}'.format(t))

        # tm.sleep(2)
        input('Press enter to continue...')

        if len(time_data) >= max_data:            
            site = int(max_data * 0.25)
            total_output -= output_avg * site
            total_state -= state_avg * site
            del time_data[0:site]
            del state_data[0:site]
            del output_data[0:site]
                        

    plt.close()
    # plt.show()

if __name__ == "__main__":    

    plant = DCMotorOpenLoop(Bm, Jm, Ki, La, Kb, Ra)

    bodePlot('DC Motor Response', plant.A(), plant.B(), plant.C(), plant.D())

    # Perturb from -100 % -> 100 %
    Bm_perturbation = Bm * (2 * random.random() - 1)
    Jm_perturbation = Jm * (2 * random.random() - 1)
    Ki_perturbation = Ki * (2 * random.random() - 1)
    La_perturbation = La * (2 * random.random() - 1)
    Kb_perturbation = Kb * (2 * random.random() - 1)
    Ra_perturbation = Ra * (2 * random.random() - 1)

    print('Perturbed by :')
    print('Bm : {} %'.format((Bm_perturbation/Bm) * 100))
    print('Jm : {} %'.format((Jm_perturbation/Jm) * 100))
    print('Ki : {} %'.format((Ki_perturbation/Ki) * 100))
    print('La : {} %'.format((La_perturbation/La) * 100))
    print('Kb : {} %'.format((Kb_perturbation/Kb) * 100))
    print('Ra : {} %'.format((Ra_perturbation/Ra) * 100))

    plant.perturbBm(Bm_perturbation)
    plant.perturbJm(Jm_perturbation)
    plant.perturbKi(Ki_perturbation)
    plant.perturbLa(La_perturbation)
    plant.perturbKb(Kb_perturbation)
    plant.perturbRa(Ra_perturbation)

    bodePlot('Perturbed', plant.A(), plant.B(), plant.C(), plant.D()) 
    # plant.unperturbAll()

    # u = 1

    # ssSimulation('Motor DC Simulation', plant.A(), plant.B(), plant.C(), plant.D(), u)

    # value = tf(1j, plant.A(), plant.B(), plant.C(), plant.D())

    # print('Value : {}'.format(value))       

    # print('A : {}'.format(A))
    # print('B : {}'.format(B))
    # print('C : {}'.format(C))
    # print('D : {}'.format(D))

    # print('A : {}'.format(A))
    # print('B1 : {}'.format(B1))
    # print('B2 : {}'.format(B2))
    # print('C1 : {}'.format(C1))
    # print('C2 : {}'.format(C2))
    # print('D11 : {}'.format(D11))
    # print('D12 : {}'.format(D12))
    # print('D21 : {}'.format(D21))
    # print('D22 : {}'.format(D22))