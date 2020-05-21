'''
    author : koseng (Lintang)
    brief : jacl Matplotlib.PyPlot Wrapper
'''

import time
import math
import numpy as np
import signal
import multiprocessing
from multiprocessing import Process, Array, Value
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

class Plotter:
    def __init__(self, n_signals, d_time, view_interval,online=True):
        self.n_signals = n_signals        
        self.max_data = int(view_interval / d_time)
        self.d_time = d_time
        self.view_interval = view_interval
        self.online = online

        self.sim_thread = Process(target = self.simulate, args = ())

        self.shared_signal_data = Array('d', np.zeros(n_signals * self.max_data))
        self.shared_time_data = Array('d', np.zeros(self.max_data))

        self.signal_data = np.frombuffer(self.shared_signal_data.get_obj()).reshape(self.max_data, n_signals)
        self.time_data = np.frombuffer(self.shared_time_data.get_obj())
        self.running = Value('b', True)
        # Default window name
        self.title = "Plotter"

        self.delay = .02

    def __del__(self):
        if(self.sim_thread.is_alive()):
            self.sim_thread.join() 

    def signalHandler(self, sig, frame):
        self.running.value = False
    
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
        print('Starting {} simulation...'.format(self.title))
        signal.signal(signal.SIGINT, self.signalHandler)        
        self.sim_thread.start()

    def pltCloseHandle(self, event):
        self.running.value = False
        print('Plotter {} closed.'.format(self.title))        

    def simulate(self):
        # Plotter setup
        plt.rcParams['savefig.facecolor'] = 'white'
        # white facecolor for temporary
        fig = plt.figure(facecolor='white', edgecolor='white')
        fig.canvas.set_window_title(self.title)
        fig.canvas.mpl_connect('close_event', self.pltCloseHandle)

        n_sp_rows, n_sp_cols = int(math.ceil(self.n_signals / 3)), min(3, self.n_signals)
        sp = []
        plot = []

        # change it to black when try to save the fig
        edge_col = 'black'
        # edge_col = 'white'
        for i in range(0, self.n_signals):
            sp.append(fig.add_subplot(n_sp_rows, n_sp_cols, i+1))
            sp[i].set_title(self.plot_name['signal{}'.format(i)], color=edge_col)            
            sp[i].grid(True)
            sp[i].set_facecolor((0.294117647, 0.294117647, 0.294117647)) # Gray
            sp[i].tick_params(color=edge_col, labelcolor=edge_col)
            sp[i].set_xlabel('Time (s)')
            sp[i].xaxis.label.set_color(edge_col)  
            for spine in sp[i].spines.values():
                spine.set_edgecolor(edge_col)

            line, = sp[i].plot(self.time_data, self.signal_data[:,i], color='r')
            plot.append(line)

        plt.tight_layout()

        # Manual adjustment
        # plt.subplots_adjust(left=0.05,bottom=0.05,right=0.99,top=0.96,wspace=0.21,hspace=0.3)
        plt.grid(True)        
        fig.canvas.draw()

        # plt.pause(0.1)
        fig.show()        
        n_data = 0

        if ~self.online:
            sp[i].set_xlim(left = self.time_data[0], right = self.time_data[-1])
            sp[i].set_ylim(bottom = min(-0.001, np.amin(self.signal_data[:,i]) * 1.75), top = max(0.001, np.amax(self.signal_data[:,i]) * 1.75))

        while self.running.value:
            n_data = n_data + 1 if n_data < self.max_data else self.max_data
            total_signal = np.array([np.sum(np.abs(self.signal_data[:,i])) for i in range(0, self.n_signals)])
            signal_avg = total_signal / len(self.signal_data) #n_data #min(n_data, self.max_data)
            
            for i in range(0, self.n_signals):
                # Replace data
                plot[i].set_ydata(self.signal_data[:,i])
                plot[i].set_xdata(self.time_data)
                # Draw into the figure
                sp[i].draw_artist(sp[i].patch)
                sp[i].draw_artist(plot[i])                
                if self.online:
                    # To maintain sight of signal in plotter
                    sp[i].set_xlim(left = max(0, self.time_data[-1] - 0.75 * self.view_interval), right = max(self.view_interval, self.time_data[-1] + 0.25 * self.view_interval))
                    # To adjust the signal interval
                    sp[i].set_ylim(bottom = min(-0.001, -signal_avg[i] * 1.75), top = max(0.001, signal_avg[i] * 1.75))
                    

            fig.canvas.draw_idle()
            fig.canvas.flush_events()

            time.sleep(self.delay)

         
