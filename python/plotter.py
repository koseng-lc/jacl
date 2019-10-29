'''
    author : koseng (Lintang)
    brief : JACL Plotter
'''

import time
import math
import numpy as np
import signal
from multiprocessing import Process, Array
import matplotlib.pyplot as plt
from scipy.linalg import expm

class Plotter:
    def __init__(self, n_signals, d_time, view_interval):
        self.n_signals = n_signals        
        self.running = True
        self.max_data = int(view_interval / d_time)
        self.d_time = d_time
        self.view_interval = view_interval

        self.sim_thread = Process(target = self.simulate, args = ())

        self.shared_signal_data = Array('d', np.zeros(n_signals * self.max_data))
        self.shared_time_data = Array('d', np.zeros(self.max_data))

        self.signal_data = np.frombuffer(self.shared_signal_data.get_obj()).reshape(self.max_data, n_signals)
        self.time_data = np.frombuffer(self.shared_time_data.get_obj())

        # Default window name
        self.title = "Plotter"

    def __del__(self):
        self.sim_thread.join()

    def signalHandler(self, sig, frame):
        self.running = False
        # self.sim_thread.join()
    
    def setData(self, _signal_data, _time_data):
        
        assert _signal_data.shape == (self.max_data, self.n_signals)
        assert _time_data.shape == (self.max_data, )

        np.copyto(self.signal_data, _signal_data)
        np.copyto(self.time_data, _time_data)

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
        # Plotter setup
        fig = plt.figure(facecolor='black', edgecolor='white')
        fig.canvas.set_window_title(self.title)
        print(self.title)
        fig.canvas.mpl_connect('close_event', self.pltCloseHandle)

        n_sp_rows, n_sp_cols = int(math.ceil(self.n_signals / 3)), min(3, self.n_signals)
        sp = []
        plot = []

        for i in range(0, self.n_signals):
            sp.append(fig.add_subplot(n_sp_rows, n_sp_cols, i+1))
            sp[i].set_title(self.plot_name['signal{}'.format(i)], color='white')
            sp[i].set_xlabel('Time (s)')
            sp[i].grid(True)
            sp[i].set_facecolor((0.294117647, 0.294117647, 0.294117647)) # Gray
            sp[i].tick_params(color='white', labelcolor='white')
            for spine in sp[i].spines.values():
                spine.set_edgecolor('white')

            line, = sp[i].plot(self.time_data, self.signal_data[:,i], color='r')
            plot.append(line)

        plt.tight_layout()

        # Manual adjustment
        # plt.subplots_adjust(left=0.05,bottom=0.05,right=0.99,top=0.96,wspace=0.21,hspace=0.3)

        plt.grid(True)        

        fig.canvas.draw()
        fig.show()

        n_data = 0

        while self.running:

            n_data += 1

            total_signal = np.array([np.sum(np.abs(self.signal_data[:,i])) for i in range(0, self.n_signals)])
            signal_avg = total_signal / min(n_data, self.max_data)

            for i in range(0, self.n_signals):

                # Replace data
                plot[i].set_ydata(self.signal_data[:,i])
                plot[i].set_xdata(self.time_data)

                # Draw into the figure
                sp[i].draw_artist(sp[i].patch)
                sp[i].draw_artist(plot[i])

                # To maintain sight of signal in plotter
                sp[i].set_xlim(left = max(0, self.time_data[-1] - 0.75 * self.view_interval), right = max(self.view_interval, self.time_data[-1] + 0.25 * self.view_interval))

                # To adjust the signal interval
                sp[i].set_ylim(bottom = min(-0.001, -signal_avg[i] * 1.75), top = max(0.001, signal_avg[i] * 1.75))

            fig.canvas.draw_idle()
            fig.canvas.flush_events()

            time.sleep(self.delay)

         
