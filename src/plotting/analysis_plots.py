from brian2 import Hz
import matplotlib.pyplot as plt
import os
import numpy as np
import scipy
import config
import params
from scipy import signal

FIGURES_DIR = config.FIGURES_DIR


def plot_auto_lfp(data):
    smoothed_data = scipy.ndimage.gaussian_filter(data, sigma=2.0)

    # plot data and smoothed data on same plot in different color
    fig, (ax1) = plt.subplots(1, 1, figsize=(10, 8), sharex=True)
    dt = params.TAU_CLOCK / params.DT_SCALING / second
    transient_samples = int(params.TRANSIENT / dt)
    data_window = data[0][transient_samples:]
    smoothed_data_window = smoothed_data[0][transient_samples:]
    ax1.plot(data_window, color='blue')
    ax1.plot(smoothed_data_window, color='orange')
    ax1.set_xlabel("Time (s)")
    plt.suptitle("x1 LFP vs Smoothed x1 (Inhibitory Population)")

    plt.savefig(os.path.join(FIGURES_DIR, "auto_lfp.png"), format="png")
    plt.show()

def plot_kop(phase_matrix):
    print(phase_matrix.shape)
    print(phase_matrix[0])
    # plot the first array in the phase matrix
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    plt.suptitle("Kop Phase For a Single Neuron")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Phase (angle in radians)")
    ax.plot(phase_matrix[0])
    plt.savefig(os.path.join(FIGURES_DIR, "kop.png"), format="png")
    plt.show()

def plot_autocorr(autocor, lag):
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    plt.suptitle("Autocorrelation")
    ax.set_xlabel("Lag (s)")
    ax.set_ylabel("Signal")
    ax.plot(lag, autocor)
    plt.savefig(os.path.join(FIGURES_DIR, "autocorr.png"), format="png")
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

def plot_power_spec(x1, x2):
    # compute mean potential
    x1_mean = np.mean(x1, axis=0)
    x2_mean = np.mean(x2, axis=0)
    x_mean = (0.8 * x1_mean) + (0.2 * x2_mean)

    fig, (ax1) = plt.subplots(1, 1, figsize=(12, 10))
    fig.suptitle("Power Spectrum (Ictal)")
    fs = 1 / (params.TAU_CLOCK/params.DT_SCALING) / Hz
    f, Pxx = signal.welch(x_mean, fs=fs)

    ax1.semilogy(f, Pxx)
    ax1.set_ylabel("Amplitude (a.u.)")
    ax1.set_xlabel("Frequency (Hz)")

    plt.savefig(os.path.join(FIGURES_DIR, "power.png"), format="png")
    plt.show()
