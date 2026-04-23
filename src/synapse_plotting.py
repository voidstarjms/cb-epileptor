from brian2 import *
from brian2tools import *
import matplotlib.pyplot as plt
import os
import numpy as np
import config
import scipy
import params as params
DATA_DIR = config.DATA_DIR
FIGURES_DIR = config.FIGURES_DIR
OUTPUT_DATA_FILE = config.OUTPUT_DATA_FILE

def plot_wpre_and_calcium(data):
    t = data['t']
    x1 = data['x1']
    wpre1 = data['syn_wpre_1']
    Ca1 = data['syn_ca_1']
    x2 = data['x2']
    wpre2 = data['syn_wpre_2']
    Ca2 = data['syn_ca_2']

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
    ax1.plot(t, wpre1[0])
    ax1.set_ylabel('Wpre')
    ax2.plot(t, Ca1[0])
    ax2.set_ylabel('Ca')
    ax3.plot(t, x1[1])
    ax3.set_ylabel('x')
    ax3.set_xlabel('time (s)')
    fig.suptitle("Synaptic dynamics with respect to x2")

    plt.savefig(os.path.join(FIGURES_DIR, "synapse_dynamics.png"), format='png')
    plt.show()