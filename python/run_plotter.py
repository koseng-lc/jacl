#!/usr/bin/env python3
import plotter

if __name__ == "__main__":
    plt = plotter.Plotter(3, 0.01, 0.1)
    plt.setPlotName({
        'signal0':'A',
        'signal1':'B',
        'signal2':'C'
    })
    plt.beginSimulation()