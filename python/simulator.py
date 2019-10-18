'''
    author : koseng (Lintang)
    brief : JACL LTI System Simulator
'''

import time
import numpy as np
import signal
import threading
import matplotlib.pyplot as plt
from scipy.linalg import expm

class Simulator:
    def __init__(self, n_states, n_inputs, n_outputs):
        self.n_states = n_states
        self.n_inputs = n_inputs
        self.n_outputs = n_outputs

        self.A = np.zeros((n_states, n_states))
        self.B = np.zeros((n_states, n_inputs))
        self.C = np.zeros((n_outputs, n_states))
        self.D = np.zeros((n_outputs, n_inputs))

        self.running = True

    def __del__(self):
        plt.close()
        # self.sim_thread.join()

    def signalHandler(self, sig, frame):
        self.running = False

    def setStateSpace(self, _A, _B, _C, _D):
        self.A = _A
        self.B = _B
        self.C = _C
        self.D = _D

    def setTitle(self, _title):
        self.title = _title

    def setDelay(self, _delay):
        self.delay = _delay

    def beginSimulation(self):
        print('Starting simulation...')
        signal.signal(signal.SIGINT, self.signalHandler)
        # self.sim_thread = threading.Thread(target=self.simulate, args=())
        # self.sim_thread.start()   
        self.simulate() 

    def simulate(self):        

        # Time interval (seconds)
        d_time = 0.0001
        view_interval = 0.01

        max_data = view_interval / d_time

        # Time
        t = 0

        # State Transition Matrix
        state_trans = expm(self.A * d_time)

        # Initial State
        init_state = np.zeros((self.n_states, 1))

        # State
        state = init_state
        last_state = state

        # Output
        output = np.zeros((self.n_outputs, 1))

        # Data
        state_data = [state[1]]
        output_data = [output[1]]
        time_data = [0]

        # Other stuff
        B_dt = self.B * d_time
        total_state = 0
        total_output = 0

        # Input
        u = np.ones((self.n_inputs, 1))
        u[0] = 0.5 # Torque
        u[1] = -12 # Voltage
        last_u = np.zeros((self.n_inputs, 1))
        
        # Plotter setup
        plt.clf()
        plt.title(self.title)

        fig = plt.figure()    

        sp1 = fig.add_subplot(211)
        sp2 = fig.add_subplot(212)

        line1, = sp1.plot(time_data, state_data, color='b')
        line2, = sp2.plot(time_data, output_data, color='r')

        fig.canvas.draw()
        fig.show()    

        sp1.set_ylim(bottom = -1.0, top = 1.0)
        sp2.set_ylim(bottom = -1.0, top = 1.0)

        while self.running:
            # Approximation of State Space solution
            term1 = np.dot(state_trans, last_state)
            term2 = np.dot(B_dt, u)
            term3 = np.dot(state_trans, np.dot(B_dt, last_u))
            term4 = np.dot(state_trans, term2 + term3) * 0.5

            state = term1 + term4

            # Output calculation based on that state approximation
            output = np.dot(self.C, state) + np.dot(self.D, u)

            last_u = u
            last_state = state      

            # Calculating state and output value
            # for calcuating their average
            total_state += np.abs(state[1])
            total_output += np.abs(output[1])

            state_data.append(state[1])
            output_data.append(output[1])
            time_data.append(t)

            # Replace data
            line1.set_ydata(state_data)
            line1.set_xdata(time_data)

            line2.set_ydata(output_data)
            line2.set_xdata(time_data)  

            sp1.draw_artist(sp1.patch)
            sp1.draw_artist(line1)

            sp2.draw_artist(sp2.patch)
            sp2.draw_artist(line2)

            # To maintain sight of signal in plotter
            sp1.set_xlim(left = max(0, t - 0.75 * view_interval), right = max(view_interval, t + 0.25 * view_interval))
            sp2.set_xlim(left = max(0, t - 0.75 * view_interval), right = max(view_interval, t + 0.25 * view_interval))    

            # To adjust the signal interval
            output_avg = total_output / len(time_data)
            state_avg = total_state / len(time_data)

            sp1.set_ylim(bottom = -state_avg - state_avg * 0.75, top = state_avg + state_avg * 0.75)
            sp2.set_ylim(bottom = -output_avg - output_avg * 0.75, top = output_avg + output_avg * 0.75)

            fig.canvas.draw()
            fig.canvas.flush_events()

            t += d_time
            
            time.sleep(self.delay)

            # Debugging part
            # print('Time Data Size : {}'.format(len(time_data)))
            # print('State Data Size : {}'.format(len(state_data)))
            # print('Output Data Size : {}'.format(len(output_data)))
            # print('State : {}'.format(state))
            # print('Output : {}'.format(output))
            # print('Time : {}'.format(t))            

            # Remove unnecessary data
            if len(time_data) >= max_data:            
                site = int(max_data * 0.25)
                total_output -= output_avg * site
                total_state -= state_avg * site
                del time_data[0:site]
                del state_data[0:site]
                del output_data[0:site]
         