'''
    author : koseng (Lintang)
    brief : JACL LTI System Simulator
'''

import time
import math
import numpy as np
import signal
from multiprocessing import Process, Array
import matplotlib.pyplot as plt
from scipy.linalg import expm

class Simulator:
    def __init__(self, n_states, n_inputs, n_outputs):
        self.n_states = n_states
        self.n_inputs = n_inputs
        self.n_outputs = n_outputs

        self.running = True
        self.d_time = 0.0001

        self.sim_thread = Process(target = self.simulate, args = ())

        self.shared_A = Array('d', np.zeros(n_states ** 2))
        self.shared_B = Array('d', np.zeros(n_states * n_inputs))
        self.shared_C = Array('d', np.zeros(n_outputs * n_states))
        self.shared_D = Array('d', np.zeros(n_outputs * n_inputs))
        self.shared_u = Array('d', np.zeros(n_inputs))
        self.shared_state_trans = Array('d', np.zeros(n_states ** 2))
        self.shared_B_dt = Array('d', np.zeros(n_states * n_inputs))

        self.A = np.frombuffer(self.shared_A.get_obj()).reshape((n_states, n_states))
        self.B = np.frombuffer(self.shared_B.get_obj()).reshape((n_states, n_inputs))
        self.C = np.frombuffer(self.shared_C.get_obj()).reshape((n_outputs, n_states))
        self.D = np.frombuffer(self.shared_D.get_obj()).reshape((n_outputs, n_inputs))
        self.u = np.frombuffer(self.shared_u.get_obj()).reshape((n_inputs, 1))
        self.state_trans = np.frombuffer(self.shared_state_trans.get_obj()).reshape((n_states, n_states))
        self.B_dt = np.frombuffer(self.shared_B_dt.get_obj()).reshape((n_states, n_inputs))

        # Default window name
        self.title = "Simulator"

    def __del__(self):
        self.sim_thread.join()

    def signalHandler(self, sig, frame):
        self.running = False
        # self.sim_thread.join()

    def setInput(self, _u):
        assert _u.shape == (self.n_inputs, 1)
        np.copyto(self.u, _u)

    def setStateSpace(self, _A, _B, _C, _D):
        assert _A.shape == (self.n_states, self.n_states)
        assert _B.shape == (self.n_states, self.n_inputs)
        assert _C.shape == (self.n_outputs, self.n_states)
        assert _D.shape == (self.n_outputs, self.n_inputs)

        np.copyto(self.A, _A)
        np.copyto(self.B, _B)
        np.copyto(self.C, _C)
        np.copyto(self.D, _D)

        # State Transition Matrix
        np.copyto(self.state_trans, expm(self.A * self.d_time))
        np.copyto(self.B_dt, self.B * self.d_time)

    def setTitle(self, _title):
        self.title = _title

    def setPlotName(self, _plot_name):
        self.plot_name = _plot_name

    def setDelay(self, _delay):
        self.delay = _delay

    def beginSimulation(self):
        print('Starting simulation...')
        signal.signal(signal.SIGINT, self.signalHandler)        
        self.sim_thread.start()

    def pltCloseHandle(self, event):
        self.running = False
        # self.sim_thread.join()
        print('Plotter closed.')

    def simulate(self):        

        # Time interval (seconds) & Maximum data        
        view_interval = 0.01
        max_data = view_interval / self.d_time

        # Time
        t = 0        

        # Initial State
        init_state = np.zeros((self.n_states, 1))

        # State
        state = init_state
        last_state = state

        # Input
        last_u = self.u

        # Output
        output = np.zeros((self.n_outputs, 1))        

        # Data

        signal_data = []

        for i in range(0, self.n_states):
            signal_data.append([state[i]])
        
        for i in range(0, self.n_inputs):
            signal_data.append([self.u[i]])

        for i in range(0, self.n_outputs):
            signal_data.append([output[i]])
        
        time_data = [0]

        # Other stuff
        total_signal = np.zeros(self.n_states + self.n_inputs + self.n_outputs)
        
        # Plotter setup
        # plt.subplots_adjust(left=.0,bottom=.0,right=.1,top=1.0)
        
        fig = plt.figure()
        fig.canvas.set_window_title(self.title)    
        fig.canvas.mpl_connect('close_event', self.pltCloseHandle)

        n_signals = self.n_states + self.n_inputs + self.n_outputs
        n_sp_rows, n_sp_cols = int(math.ceil(n_signals / 3)), min(3, n_signals)
        sp = []
        plot = []

        for i in range(0, n_signals):
            sp.append(fig.add_subplot(n_sp_rows, n_sp_cols, i+1))
            sp[i].set_title(self.plot_name['signal{}'.format(i)])
            sp[i].set_xlabel('Time (s)')
            sp[i].grid(True)

            color = 'b'
            if i < self.n_states:
                color = 'r'
            elif i < self.n_states + self.n_inputs:
                color = 'g'

            line, = sp[i].plot(time_data, signal_data[i], color=color)
            plot.append(line)

        plt.tight_layout()
        plt.grid(True)

        fig.canvas.draw()
        fig.show()

        while self.running:
            # Approximation of State Space solution
            term1 = np.dot(self.state_trans, last_state)
            term2 = np.dot(self.B_dt, self.u)
            term3 = np.dot(self.state_trans, np.dot(self.B_dt, last_u))
            term4 = np.dot(self.state_trans, term2 + term3) * 0.5

            state = term1 + term4

            # Output calculation based on that state approximation
            output = np.dot(self.C, state) + np.dot(self.D, self.u)

            last_u = self.u
            last_state = state                  

            pres_signal = np.concatenate((state, np.concatenate((self.u, output), axis = 0)), axis = 0)

            time_data.append(t)

            time_data_len = len(time_data)

            signal_avg = []

            for i in range(0, n_signals):

                # Summing signal magnitude
                # for calcuating their average
                total_signal[i] += np.abs(pres_signal[i])

                # Replace data
                signal_data[i].append(pres_signal[i])
                plot[i].set_ydata(signal_data[i])
                plot[i].set_xdata(time_data)

                # Draw into the figure
                sp[i].draw_artist(sp[i].patch)
                sp[i].draw_artist(plot[i])

                # To maintain sight of signal in plotter
                sp[i].set_xlim(left = max(0, t - 0.75 * view_interval), right = max(view_interval, t + 0.25 * view_interval))

                # To adjust the signal interval
                signal_avg.append(total_signal[i] / time_data_len)
                sp[i].set_ylim(bottom = min(-0.001, -signal_avg[i] * 1.75), top = max(0.001, signal_avg[i] * 1.75))

            fig.canvas.draw_idle()
            fig.canvas.flush_events()

            t += self.d_time
            
            time.sleep(self.delay)

            # Debugging part
            # print('Time Data Size : {}'.format(len(time_data)))
            # print('State Data Size : {}'.format(len(state_data)))
            # print('Output Data Size : {}'.format(len(output_data)))
            # print('State : {}'.format(state))
            # print('Output : {}'.format(output))
            # print('Time : {}'.format(t))            

            # Remove unnecessary old data
            if len(time_data) >= max_data:

                site = int(max_data * 0.25)

                for i in range(0, n_signals):

                    total_signal[i] -= signal_avg[i] * site
                    del signal_data[i][0:site]

                del time_data[0:site]
         