from brian2 import *
from brian2tools import *
import matplotlib.pyplot as plt
import os
import numpy as np
import config

DATA_DIR = config.DATA_DIR
FIGURES_DIR = config.FIGURES_DIR
OUTPUT_DATA_FILE = config.OUTPUT_DATA_FILE


# kwargs collects args into a dict, allows flexible arguments to be passed in
def save_data(filename, **kwargs):
    if not os.path.exists(DATA_DIR):
        os.makedirs(DATA_DIR)
    
    np.savez(os.path.join(DATA_DIR, filename), **kwargs)

def standard_plot():
    pass

def raster(fig_name: str, population_name: str, t, x, spike_matrix, num_cells, sim_duration):
    
    x_mean = np.mean(x, axis=0)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(30, 10), sharex=True)
    fig.suptitle(f'All {population_name} Variables - All Neurons Averaged')

    # Plot the averaged data instead of just neuron 0
    ax1.plot(t, x_mean)
    ax1.set_ylabel("Mean x") 
    
    # configure main raster plot
    raster = ax2.imshow(spike_matrix, interpolation='none', aspect='auto',
                   origin='lower', extent=[0, sim_duration, 0, num_cells])

    ax2.set_xlabel('Time (s)', fontsize=12)
    ax2.set_ylabel('Neuron index', fontsize=12)
    ax2.set_title(f'{population_name} Spike Raster (Spike Count)', fontsize=14)

    # config colorbar
    cbar = fig.colorbar(raster, ax=ax2, location='bottom', aspect=50)
    cbar.minorticks_on()

    # save plot
    plt.savefig(os.path.join(FIGURES_DIR, f"{fig_name}_raster.png"), format='png')
    plt.show()



def plot_hr_single(t, x1, y1, z1, I_syn_inter_1):
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(30, 10), sharex=True)
    fig.suptitle("All Hindmarsh Rose Variables - One Neuron")
    ax1.plot(t, x1[0])
    ax1.set_ylabel("Neuron 0 x")
    ax2.plot(t, y1[0])
    ax2.set_ylabel("Neuron 0 y")
    ax3.plot(t, z1[0])
    ax3.set_ylabel("Neuron 0 z")
    ax4.plot(t, I_syn_inter_1[0])
    ax4.set_ylabel("Neuron 0 I_{syn, inter}")
    ax4.set_xlabel("Time (s)")
    plt.savefig(os.path.join(FIGURES_DIR, "pop1_single_neuron.png"), format="png")
    plt.show()


def plot_hr_mean(t, x1, y1, z1):
    # Calculate the mean across all neurons (axis=0)
    x1_mean = np.mean(x1, axis=0)
    y1_mean = np.mean(y1, axis=0)
    z1_mean = np.mean(z1, axis=0)

    # All pop 1 variables (Figure 2 - Now Averaged)
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(30, 10), sharex=True)
    fig.suptitle("All Hindmarsh Rose Variables - All Neurons Averaged")

    # Plot the averaged data instead of just neuron 0
    ax1.plot(t, x1_mean)
    ax1.set_ylabel("Mean x1") # Updated label
    ax2.plot(t, y1_mean)
    ax2.set_ylabel("Mean y1") # Updated label
    ax3.plot(t, z1_mean)
    ax3.set_ylabel("Mean z1") # Updated label

    ax3.set_xlabel("Time (s)")

    plt.savefig(os.path.join(FIGURES_DIR, "pop1_mean_neurons.png"), format="png")
    plt.show()

def plot_ml_single(t, x2, n2):
    # All pop2 variables
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(30, 10), sharex=True)
    fig.suptitle("All Morris Lecar Variables - One Neurons")
    ax1.plot(t, x2[0])
    ax1.set_ylabel("Neuron 0 x")
    ax2.plot(t, n2[0])
    ax2.set_ylabel("Neuron 0 n")
    ax2.set_xlabel("Time (s)")
    plt.savefig(os.path.join(FIGURES_DIR, "pop2_single_neuron.png"), format="png")
    plt.show()

def plot_ml_mean():
    pass



def plot_both(t, x1, x2):
    # One neuron from both pops 
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.suptitle("One Neuron From Each Population")
    ax1.plot(t, x1[0])
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("x1")
    ax2.plot(t, x2[0])
    ax2.set_xlabel("Time (s)")
    ax2.set_ylabel("x2")


    #plt.savefig("figures/interictal_pop2_r4e-5_10s.png", format="png")
    plt.savefig(os.path.join(FIGURES_DIR, "pop1_and_pop2_single.png"), format="png")
    plt.show()


def plot_both_avg(t, x1, y1, z1, x2, n):
    # neurons averaged from both pops 
    x1_mean = np.mean(x1, axis=0)
    y1_mean = np.mean(y1, axis=0)
    z1_mean = np.mean(z1, axis=0)
    x2_mean = np.mean(x2, axis=0)
    n_mean = np.mean(n, axis=0)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.suptitle("Mean Neuron From Each Population")
    ax1.plot(t, x1_mean)
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("x1")
    ax2.plot(t, x2_mean)
    ax2.set_xlabel("Time (s)")
    ax2.set_ylabel("x2")


    #plt.savefig("figures/interictal_pop2_r4e-5_10s.png", format="png")
    plt.savefig(os.path.join(FIGURES_DIR, "pop1_and_pop2_single.png"), format="png")
    plt.show()
    
def plot_mean_potential():
    pass
    # pop1_mean = np.mean(x1, axis=0)
    # pop2_mean = np.mean(x2, axis=0)
    # mean_potential = 0.8 * pop1_mean + 0.2 * pop2_mean
    # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 8))
    # ax1.plot(t, pop2_mean)
    # ax1.set_xlabel("Time (s)")
    # ax1.set_ylabel("Weighted mean potential (a.u.)")

def plot_power_spec():
    pass
    
    # fs = 1 / defaultclock.dt / Hz
    # f, Pxx = signal.welch(mean_potential, fs=fs)
    # ax2.semilogy(f, Pxx)
    # ax2.set_xlabel("Frequency (Hz)")
    # ax2.set_ylabel("Amplitude (a.u.)")
    # plt.savefig("figures/interictal_power_spectrum_hi_r.png", format="png")
    # plt.show()

    # fig = plt.figure()
    # f, ts, Sxx = scipy.signal.spectrogram(mean_potential, fs)
    # fig = plt.pcolormesh(ts, f, Sxx, shading='gouraud')
