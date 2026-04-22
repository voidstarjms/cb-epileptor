"""Signal-analysis plots: smoothed LFP, KOP phase, autocorrelation, power spectrum."""
from brian2 import *
import matplotlib.pyplot as plt
import os
import numpy as np
import scipy
from scipy import signal
from typing import Dict, Any


def plot_auto_lfp(filepaths: Any, params_dict: Dict, data) -> None:
    """Overlay raw and Gaussian-smoothed traces for neuron 0, past the transient.
        INPUT:
            filepaths: FilePaths with figures_dir.
            params_dict: needs TAU_CLOCK, DT_SCALING, TRANSIENT.
            data: (num_neurons, time) signal.
    """
    smoothed_data = scipy.ndimage.gaussian_filter(data, sigma=2.0)

    # plot data and smoothed data on same plot in different color
    fig, (ax1) = plt.subplots(1, 1, figsize=(10, 8), sharex=True)
    dt = params_dict['TAU_CLOCK'] / params_dict['DT_SCALING'] / second
    transient_samples = int(params_dict['TRANSIENT'] / dt)
    data_window = data[0][transient_samples:]
    smoothed_data_window = smoothed_data[0][transient_samples:]
    ax1.plot(data_window, color='blue')
    ax1.plot(smoothed_data_window, color='orange')
    ax1.set_xlabel("Time (s)")
    plt.suptitle("x1 LFP vs Smoothed x1 (Inhibitory Population)")

    plt.savefig(os.path.join(filepaths.figures_dir, "auto_lfp.png"), format="png")
    plt.show()

def plot_kop(filepaths: Any, phase_matrix: np.ndarray) -> None:
    """Plot the interpolated phase trace of the first neuron.
        INPUT:
            filepaths: FilePaths with figures_dir.
            phase_matrix: (num_neurons, time) phase in radians, from synch.compute_phase.
    """
    # plot the first array in the phase matrix
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    plt.suptitle("Kop Phase For a Single Neuron")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Phase (angle in radians)")
    ax.plot(phase_matrix[0])
    plt.savefig(os.path.join(filepaths.figures_dir, "kop.png"), format="png")
    plt.show()

def plot_autocorr(filepaths: Any, autocorr: np.ndarray, lag: np.ndarray) -> None:
    """Plot autocorrelation vs. lag.
        INPUT:
            filepaths: FilePaths with figures_dir.
            autocorr: autocorrelation values.
            lag: lag values, same length as autocorr.
    """
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    plt.suptitle("Autocorrelation")
    ax.set_xlabel("Lag (s)")
    ax.set_ylabel("Signal")
    ax.plot(lag, autocorr)
    plt.savefig(os.path.join(filepaths.figures_dir, "autocorr.png"), format="png")
    plt.show()


def plot_mean_potential():
    """Not implemented."""
    pass
    # pop1_mean = np.mean(x1, axis=0)
    # pop2_mean = np.mean(x2, axis=0)
    # mean_potential = 0.8 * pop1_mean + 0.2 * pop2_mean
    # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 8))
    # ax1.plot(t, pop2_mean)
    # ax1.set_xlabel("Time (s)")
    # ax1.set_ylabel("Weighted mean potential (a.u.)")

def plot_power_spec(filepaths: Any, params_dict: Dict, x1: np.ndarray, x2: np.ndarray) -> None:
    """Welch power spectrum of the 80/20 HR/ML weighted-mean signal.
        INPUT:
            filepaths: FilePaths with figures_dir.
            params_dict: needs TAU_CLOCK and DT_SCALING to set the sampling rate.
            x1: (num_cells, time) HR x.
            x2: (num_cells, time) ML x.
    """
    # compute mean potential
    x1_mean = np.mean(x1, axis=0)
    x2_mean = np.mean(x2, axis=0)
    x_mean = (0.8 * x1_mean) + (0.2 * x2_mean)

    fig, (ax1) = plt.subplots(1, 1, figsize=(12, 10))
    fig.suptitle("Power Spectrum (Ictal)")
    fs = 1 / (params_dict['TAU_CLOCK']/params_dict['DT_SCALING']) / Hz
    f, Pxx = signal.welch(x_mean, fs=fs)

    ax1.semilogy(f, Pxx)
    ax1.set_ylabel("Amplitude (a.u.)")
    ax1.set_xlabel("Frequency (Hz)")

    plt.savefig(os.path.join(filepaths.figures_dir, "power.png"), format="png")
    plt.show()
