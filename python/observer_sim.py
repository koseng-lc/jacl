'''
    author : koseng (Lintang)
    brief : Full-Order Observer Simulator
'''

import time
import math
import numpy as np
import signal
from multiprocessing import Process, Array
import matplotlib.pyplot as plt
from scipy.linalg import expm

class ObserverSimulator:
    def __init__(self, n_states, n_inputs, n_outputs):
        self.n_states = n_states
        self.n_inputs = n_inputs
        self.n_outputs = n_outputs

        self.running = True
        self.d_time = 0.0001

        self.sim_thread = Process(target = self.simulate, args = ())

        self.shared_A_hat = Array('d', np.zeros(n_states ** 2))
        self.shared_K = Array('d', np.zeros(n_states * n_outputs))
        self.shared_H = Array('d', np.zeros(n_states * n_inputs))
        self.shared_C = Array('d', np.zeros(n_outputs * n_states))
        self.shared_D = Array('d', np.zeros(n_outputs * n_inputs))
        self.shared_u = Array('d', np.zeros(n_inputs))
        self.shared_state_trans = Array('d', np.zeros(n_states ** 2))
        
        # Shared (State Transition Matrix + I (Identity)) * d_time
        self.shared_st_I_dt = Array('d', np.zeros(n_states ** 2))

        self.A_hat = np.frombuffer(self.shared_A_hat.get_obj()).reshape((n_states, n_states))
        self.K = np.frombuffer(self.shared_K.get_obj()).reshape((n_states, n_outputs))
        self.H = np.frombuffer(self.shared_H.get_obj()).reshape((n_states, n_inputs))
        self.C = np.frombuffer(self.shared_C.get_obj()).reshape((n_outputs, n_states))
        self.D = np.frombuffer(self.shared_D.get_obj()).reshape((n_outputs, n_inputs))
        self.u = np.frombuffer(self.shared_u.get_obj()).reshape((n_inputs, 1))
        self.state_trans = np.frombuffer(self.shared_state_trans.get_obj()).reshape((n_states, n_states))
        self.st_I_dt = np.frombuffer(self.shared_st_I_dt.get_obj()).reshape((n_states, n_states))

        # Default window name
        self.title = "Observer Simulator"

    def __del__(self):
        self.sim_thread.join()

    def signalHandler(self, sig, frame):
        self.running = False
        # self.sim_thread.join()

    def setInput(self, _u):
        assert _u.shape == (self.n_inputs, 1)
        np.copyto(self.u, _u)

    def setStateSpace(self, _A_hat, _K, _H, _C, _D):
        assert _A_hat.shape == (self.n_states, self.n_states)
        assert _K.shape == (self.n_states, self.n_outputs)
        assert _H.shape == (self.n_states, self.n_inputs)
        assert _C.shape == (self.n_outputs, self.n_states)
        assert _D.shape == (self.n_outputs, self.n_inputs)

        np.copyto(self.A_hat, _A_hat)
        np.copyto(self.K, _K)
        np.copyto(self.H, _H)
        np.copyto(self.C, _C)
        np.copyto(self.D, _D)

        # State Transition Matrix
        np.copyto(self.state_trans, expm(self.A_hat * self.d_time))
        np.copyto(self.st_I_dt, (np.eye(self.n_states) + self.state_trans) * self.d_time)

    def setTitle(self, _title):
        self.title = _title

    def setPlotName(self, _plot_name):
        self.plot_name = _plot_name

    def setDelay(self, _delay):
        self.delay = _delay

    def beginSimulation(self):
        print('Starting Observer simulation...')
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
        fig = plt.figure(facecolor='black', edgecolor='white')
        fig.canvas.set_window_title(self.title)    
        fig.canvas.mpl_connect('close_event', self.pltCloseHandle)

        n_signals = self.n_states + self.n_inputs + self.n_outputs
        n_sp_rows, n_sp_cols = int(math.ceil(n_signals / 3)), min(3, n_signals)
        sp = []
        plot = []

        for i in range(0, n_signals):
            sp.append(fig.add_subplot(n_sp_rows, n_sp_cols, i+1))
            sp[i].set_title(self.plot_name['signal{}'.format(i)], color='white')
            sp[i].set_xlabel('Time (s)')
            sp[i].grid(True)
            sp[i].set_facecolor((0.294117647, 0.294117647, 0.294117647)) # Gray
            sp[i].tick_params(color='white', labelcolor='white')
            for spine in sp[i].spines.values():
                spine.set_edgecolor('white')

            color = 'b'
            if i < self.n_states:
                color = 'r'
            elif i < self.n_states + self.n_inputs:
                color = 'g'

            line, = sp[i].plot(time_data, signal_data[i], color=color)
            plot.append(line)

        plt.tight_layout()

        # Manual adjustment
        plt.subplots_adjust(left=0.05,bottom=0.05,right=0.99,top=0.96,wspace=0.21,hspace=0.3)

        plt.grid(True)

        fig.canvas.draw()
        fig.show()

        while self.running:
            
            # Approximation of State Space solution
            term1 = np.dot(self.state_trans, last_state)
            term2 = np.dot(self.K, output)
            term3 = np.dot(self.H, self.u)
            term4 = np.dot(self.st_I_dt, term2 + term3) * .5

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
         