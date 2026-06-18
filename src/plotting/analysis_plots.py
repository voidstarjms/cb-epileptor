"""Signal-analysis plots: smoothed LFP, KOP phase, autocorrelation, power spectrum.

Each plot writes a PNG under filepaths.figures_dir.
"""
from brian2 import *
import matplotlib.pyplot as plt
import os
import numpy as np
import scipy
from scipy import signal
from typing import Dict, Any

import plotting.style  # noqa: F401


def plot_auto_lfp(filepaths: Any, params_dict: Dict, data) -> None:
    """Overlay raw and Gaussian-smoothed traces for neuron 0, past the transient.

    Args:
        filepaths (Any): FilePaths with figures_dir.
        params_dict (Dict): Needs TAU_CLOCK, DT_SCALING, TRANSIENT.
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
    plt.suptitle("x LFP vs Smoothed x (Inhibitory Population)")

    plt.savefig(os.path.join(filepaths.figures_dir, "auto_lfp.png"), format="png")
    plt.show()

def plot_kop(filepaths: Any, phase_matrix: np.ndarray) -> None:
    """Plot the interpolated phase trace of the first neuron.

    Args:
        filepaths (Any): FilePaths with figures_dir.
        phase_matrix (np.ndarray): (num_neurons, time) phase in radians, from
            synch.compute_phase.
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

    Args:
        filepaths (Any): FilePaths with figures_dir.
        autocorr (np.ndarray): Autocorrelation values.
        lag (np.ndarray): Lag values, same length as autocorr.
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

def plot_power_spec(filepaths: Any, params_dict: Dict, x1: np.ndarray, x2: np.ndarray, fmax: float = 100.0) -> None:
    """Welch power spectrum of the 80/20 HR/ML weighted-mean signal.

    Args:
        filepaths (Any): FilePaths with figures_dir.
        params_dict (Dict): Needs TAU_CLOCK and DT_SCALING to set the sampling rate.
        x1 (np.ndarray): (num_cells, time) HR x.
        x2 (np.ndarray): (num_cells, time) ML x.
        fmax (float): Upper frequency limit for display in Hz. Default 100 Hz.
    """
    x1_mean = np.mean(x1, axis=0)
    x2_mean = np.mean(x2, axis=0)
    x_mean = (0.8 * x1_mean) + (0.2 * x2_mean)

    fs = float(1 / (params_dict['TAU_CLOCK'] / params_dict['DT_SCALING']) / Hz)
    f, Pxx = signal.welch(x_mean, fs=fs)

    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    fig.suptitle("Power Spectrum")
    ax.semilogy(f, Pxx)
    ax.set_xlim(0, fmax)
    ax.set_ylabel("Power (a.u.)")
    ax.set_xlabel("Frequency (Hz)")

    plt.savefig(os.path.join(filepaths.figures_dir, "power.png"), format="png")
    plt.show()


def plot_spectrogram(filepaths: Any, params_dict: Dict, x1: np.ndarray, x2: np.ndarray,
                     t: np.ndarray = None, x0_t: np.ndarray = None, ce_t: np.ndarray = None,
                     fmax: float = 100.0) -> None:
    """Time-frequency spectrogram of the 80/20 HR/ML weighted-mean LFP.

    Optionally adds an x0/CE schedule panel below when t, x0_t, and ce_t
    are provided.

    Args:
        filepaths (Any): FilePaths with figures_dir.
        params_dict (Dict): Needs TAU_CLOCK, DT_SCALING, and TRANSIENT.
        x1 (np.ndarray): (num_cells, time) HR x.
        x2 (np.ndarray): (num_cells, time) ML x.
        t (np.ndarray): Time vector matching x1/x2 (seconds). Required for param panel.
        x0_t (np.ndarray): Recorded x0 schedule over time. If None, param panel is skipped.
        ce_t (np.ndarray): Recorded CE schedule over time. If None, param panel is skipped.
        fmax (float): Upper frequency limit for display in Hz. Default 100 Hz.
    """
    x_lfp = 0.8 * np.mean(x1, axis=0) + 0.2 * np.mean(x2, axis=0)

    fs = float(1 / (params_dict['TAU_CLOCK'] / params_dict['DT_SCALING']) / Hz)

    has_params = t is not None and x0_t is not None and ce_t is not None
    n_axes = 2 if has_params else 1
    height_ratios = [4, 1] if has_params else [1]

    fig, axes = plt.subplots(n_axes, 1, figsize=(20, 6 + 2 * has_params), sharex=True,
                             layout='constrained',
                             gridspec_kw={'height_ratios': height_ratios})
    ax_spec = axes[0] if has_params else axes
    fig.suptitle("LFP Spectrogram")
    fig.set_constrained_layout_pads(w_pad=0.1, h_pad=0.1, wspace=0.02, hspace=0.02)

    _, _, _, im = ax_spec.specgram(x_lfp, Fs=fs, NFFT=int(fs), noverlap=int(fs) // 2,
                                   xextent=(params_dict['TRANSIENT'],
                                            params_dict['TRANSIENT'] + len(x_lfp) / fs),
                                   cmap='viridis')
    ax_spec.set_ylim(0, fmax)
    ax_spec.set_ylabel("Frequency (Hz)")
    fig.colorbar(im, ax=ax_spec, location='right', aspect=25, pad=0.001, label="Power (dB)")

    if has_params:
        ax_p = axes[1]
        ax_p.plot(t, x0_t, label='x0', color='blue')
        ax_p.plot(t, ce_t, label='CE', color='orange')
        ax_p.set_ylabel("Value")
        ax_p.legend(loc='upper right', fontsize=10)
        ax_p.set_xlabel("Time (s)")
    else:
        ax_spec.set_xlabel("Time (s)")

    plt.savefig(os.path.join(filepaths.figures_dir, "spectrogram.png"), format="png")
    plt.show()
