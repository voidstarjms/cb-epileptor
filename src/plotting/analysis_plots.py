"""Signal-analysis plots: smoothed LFP, KOP phase, autocorrelation, power spectrum.

Each plot writes a PNG under filepaths.figures_dir.
"""
from brian2 import *
import matplotlib.pyplot as plt
import os
import numpy as np
import scipy
from scipy import signal, stats
from typing import Dict, Any
from matplotlib.backends.backend_pdf import PdfPages

import plotting.style  # noqa: F401

EPHYS_FS = 10000
EPHYS_POWER_NORMALIZE = 10e10

def plot_auto_lfp(filepaths: Any, params_dict: Dict, data) -> None:
    """Overlay raw and Gaussian-smoothed traces for neuron 0, past the transient.

    Args:
        filepaths (Any): FilePaths with figures_dir.
        params_dict (Dict): Needs TAU_CLOCK, DT_SCALING, TRANSIENT.
        data: (num_neurons, time) signal.
    
    Returns: void
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

def plot_kop(filepaths: Any, phase_matrix: np.ndarray, t_start: float, t_end: float) -> None:
    """Plot the interpolated phase trace of the first neuron.

    Args:
        filepaths (Any): FilePaths with figures_dir.
        phase_matrix (np.ndarray): (num_neurons, time) phase in radians, from
            synch.compute_phase.
        t_start (float): Start time of the phase grid in seconds (TRANSIENT).
        t_end (float): End time of the phase grid in seconds (TRANSIENT + SIM_DURATION).

    Returns: void
    """
    # compute_phase's grid includes both endpoints, so linspace reproduces it
    t = np.linspace(t_start, t_end, phase_matrix.shape[1])

    # plot the first array in the phase matrix
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    plt.suptitle("Kop Phase For a Single Neuron")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Phase (angle in radians)")
    ax.plot(t, phase_matrix[0])
    ax.set_xlim(t_start, t_end)
    plt.savefig(os.path.join(filepaths.figures_dir, "kop.png"), format="png")
    plt.show()

def plot_autocorr(filepaths: Any, autocorr: np.ndarray, lag: np.ndarray) -> None:
    """Plot autocorrelation vs. lag.

    Args:
        filepaths (Any): FilePaths with figures_dir.
        autocorr (np.ndarray): Autocorrelation values.
        lag (np.ndarray): Lag values, same length as autocorr.
    
    Returns: void
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

def _make_power_spec_plot(lfp_list : list, fs : float, fmax : float,
                          bounds : tuple = (None, None),
                          labels : list = None):
    """Private core function for generating power spectral density plots.
    
    Args:
        fig_dir (str): Output directory for figure.
        lfp_list (list[np.ndarray]): list of LFP arrays on which to compute PSD.
        fs (float): Sample frequency of LFPs.
        fmax (float): Max frequency to plot.
        fname (str): Output file name.
        bounds (tuple[None | float, None | float]):
            min_bound: lower bound of y-axis (dB).
            max_bound: upper bound of y-axis (dB).
        labels (list): L

    Returns: void
    """
    
    if labels != None and len(labels) != len(lfp_list):
        print(f"_make_power_spec_plot: Label count does not match LFP count. Label count is {len(labels)} but LFP count is {len(lfp_list)}")
        return
    
    f_list = []
    Pxx_list = []
    for i in lfp_list:
        f, Pxx = signal.welch(i, fs=fs, nperseg=16384, scaling='spectrum')
        f_list.append(f)
        Pxx_list.append(Pxx)

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    if labels == None:
        for f, Pxx in zip(f_list, Pxx_list):
            ax.semilogy(f, Pxx)
    else:
        for f, Pxx, l in zip(f_list, Pxx_list, labels):
            ax.semilogy(f, Pxx, label=l)
    
    ax.set_xlim(0, fmax)
    min_bound, max_bound = bounds
    if min_bound != None:
        ax.set_ylim(bottom=min_bound)
    if max_bound != None:
        ax.set_ylim(top=max_bound)

    ax.set_ylabel("Power (dB)")
    ax.set_xlabel("Frequency (Hz)")
    plt.xticks(np.arange(0, fmax, fmax / 10))
    if labels != None:
        plt.legend()

    return fig

def _output_power_spec_plot(fig_dir : str, lfp_list : list, fs : float, fmax : float,
                          fname : str = "power", bounds : tuple = (None, None),
                          labels : list = None):
    fig = plt.figure(_make_power_spec_plot(lfp_list, fs, fmax,
                                     bounds, labels))
    fig.suptitle("Power Spectrum")
    
    plt.savefig(os.path.join(fig_dir, fname + ".png"), format="png")
    plt.show()

def plot_power_spec(filepaths: Any, x1: np.ndarray, x2: np.ndarray, fs: float,
                    fmax: float = 40.0, fname : str = "power") -> None:
    """Welch power spectrum of the 80/20 EPOP/IPOP weighted-mean signal.

    Args:
        filepaths (Any): FilePaths with figures_dir.
        x1 (np.ndarray): (num_cells, time) EPOP x.
        x2 (np.ndarray): (num_cells, time) IPOP x.
        fs (float): Sampling rate of x1/x2 in Hz (see synch_metrics.sample_rate).
        fmax (float): Upper frequency limit for display in Hz. Default 100 Hz.

    Returns: void
    """
    x1_mean = np.mean(x1, axis=0)
    x2_mean = np.mean(x2, axis=0)
    x_mean = (0.8 * x1_mean) + (0.2 * x2_mean)

    _output_power_spec_plot(filepaths.figures_dir, [x_mean], fs, fmax, fname=fname)

def plot_ephys_power_spec(fig_dir : str, lfp_list : list, fmax : float = 100.0,
                          fname : str = "ephys_power", labels : list = None,
                          bounds=(None, None)):
    _output_power_spec_plot(fig_dir, lfp_list, EPHYS_FS, fmax, fname=fname,
                          bounds=bounds, labels=labels)

def ephys_power_spec_pdf(lfp_list : list, fmax : float = 100.0,
                        outfile : str = "ephys_mean_power", labels : list = None):
    with PdfPages(outfile) as pdf:
        for i, s in enumerate(lfp_list):
            if s is not np.nan:
                fig = plt.figure(_make_power_spec_plot(lfp_list[i:i+1], EPHYS_FS, fmax,
                                                 #bounds=(10e-17, 10e-11),
                                                 labels=labels))
                fig.suptitle(f"Power Spectrum (Sweep {i+1})")

                pdf.savefig()
                plt.close()

def _normalize_power_spec(freq : list, spec_list : list):
    freq_count = len(freq)
    noise_floor_start_idx = 0
    while noise_floor_start_idx < freq_count and freq[noise_floor_start_idx] < 70:
        noise_floor_start_idx += 1
    noise_floor_end_idx = noise_floor_start_idx
    while noise_floor_end_idx < freq_count and freq[noise_floor_end_idx] < 100:
        noise_floor_end_idx += 1

    power_accum = 0
    nan_sweep_count = 0
    for s in spec_list:
        if s is not np.nan:
            power_accum += np.sum(np.square(s))
            noise_floor = (np.mean(
                s[noise_floor_start_idx:noise_floor_end_idx]))
            s -= noise_floor
        else:
            nan_sweep_count += 1
    spec_count = len(spec_list)
    mean_total_power = np.sqrt(power_accum / (spec_count - nan_sweep_count))

    for s in spec_list:
        if s is not np.nan:
            s /= mean_total_power

    return spec_list

def _compute_specgram_list(lfp_list : list, fs=EPHYS_FS):
    spec_list = []
    freq = None
    for s in lfp_list:
        if s is not np.nan:
            freq, Pxx = signal.welch(s, fs=fs, nperseg=16384, scaling='spectrum')
            spec_list.append(np.log10(Pxx * EPHYS_POWER_NORMALIZE))
        else:
            spec_list.append(np.nan)

    return spec_list, freq

def plot_ephys_mean_power_spec(fig_dir : str, lfp_list : list, prepep_sweep_counts : list[int],
                                fmax : float = 100.0, fname : str = "ephys_mean_power",
                                labels : list = None):
    spec_list_list = []
    freqs = None
    for lfp in lfp_list:
        spec_list, freqs = _compute_specgram_list(lfp)
        #spec_list_list.append(spec_list)
        spec_list_list.append(_normalize_power_spec(freqs, spec_list))
    min_prepep_sweep_count = np.min(prepep_sweep_counts)
    min_postpep_sweep_count = np.min([len(spec_list_list[i]) - prepep_sweep_counts[i]
                                      for i in range(len(spec_list_list))])
    for i, s in enumerate(spec_list_list):
        spec_list_list[i] = s[prepep_sweep_counts[i] - min_prepep_sweep_count:
              prepep_sweep_counts[i] + min_postpep_sweep_count]
    
    spec_list_array = np.asarray(spec_list_list)
    spec_list_array = np.swapaxes(spec_list_array, 0, 1)
    outfile = os.path.join(fig_dir, fname)

    # TODO factor this out, the same code is used elsewhere
    n_expts = len(spec_list_list)
    with PdfPages(outfile) as pdf:
        for i, s in enumerate(list(spec_list_array)):
            if s is not np.nan:
                fig = plt.figure()
                ax = plt.gca()
                mu = s.mean(axis=0)
                ax.plot(freqs, mu)

                # Compute and plot confidence interval
                sem = s.std(axis=0, ddof=1) / np.sqrt(n_expts)
                half_width = stats.t.ppf(0.975, df=n_expts - 1) * sem
                ax.fill_between(freqs, mu - half_width, mu + half_width,
                            color="C0", alpha=0.25, linewidth=0, label="95% CI")

                plt.suptitle("PEP"+("+" if i >= min_prepep_sweep_count else "")
                             +f"{i-min_prepep_sweep_count}", fontsize=25)
                ax.set_xlim(0, fmax)
                ax.set_ylim(-0.5e-3, 2e-3)
            
                ax.set_ylabel("Log-Power (dB)", fontsize=15)
                ax.set_xlabel("Frequency (Hz)", fontsize=15)
                ax.xaxis.set_tick_params('both', labelsize=10)
                ax.yaxis.set_tick_params('both', labelsize=10)
                plt.xticks(np.arange(0, fmax, fmax / 10))

                plt.tight_layout()
                pdf.savefig()
                plt.close()

def ephys_normalized_power_pdf(lfp_list : list, fmax : float = 100.0,
                               outfile : str = "ephys_power_normalized"):
    spec_list, freq = _compute_specgram_list(lfp_list)    
    spec_list = _normalize_power_spec(freq, spec_list)

    with PdfPages(outfile) as pdf:
        for i, s in enumerate(spec_list):
            if s is not np.nan:
                fig = plt.figure(figsize=(12, 8))
                ax = plt.gca()
                
                ax.plot(freq, s)

                fig.suptitle(f"Sweep {i+1}")
                ax.set_xlim(0, fmax)
                ax.set_ylim(-1e-3, 5e-3)
            
                ax.set_ylabel("Log-Power (dB)")
                ax.set_xlabel("Frequency (Hz)")
                plt.xticks(np.arange(0, fmax, fmax / 10))

                pdf.savefig()
                plt.close()

def ephys_normalized_power(lfp_list : list, fmax : float = 100.0,
                            outfile : str = "ephys_power_normalized",
                            sweep : int = -1):
    # Ignore all sweeps except the requested one if sweep != -1
    if sweep != -1:
        lfp_list = lfp_list[sweep]
    
    spec_list, freq = _compute_specgram_list(lfp_list)
    spec_list = _normalize_power_spec(freq, spec_list)
    spec_count = len(spec_list)

    fig = plt.figure(figsize=(10, 8))
    ax = plt.gca()
    cmap = plt.cm.plasma
    col_list = cmap(np.linspace(0, 1, spec_count))
    norm = matplotlib.colors.Normalize(vmin=0, vmax=spec_count)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    plot_alpha = 0.4 if sweep == -1 else 1
    for i, s in enumerate(spec_list):
        if s is not np.nan:
            ax.plot(freq, s, c=col_list[i], alpha=plot_alpha)
    fig.colorbar(sm, ax=ax, location='right', aspect=25, pad=0.001, label="Sweep number")

    fig.suptitle(f"Power Spectrum")
    ax.set_xlim(0, fmax)
    ax.set_ylim(-1e-3, 5e-3)
    ax.set_facecolor("#000000")

    ax.set_ylabel("Log-Power (dB)")
    ax.set_xlabel("Frequency (Hz)")
    plt.xticks(np.arange(0, fmax, fmax / 10))

    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()

def plot_spectrogram(filepaths: Any, params_dict: Dict, x1: np.ndarray, x2: np.ndarray,
                     fs: float, t: np.ndarray = None, x0_t: np.ndarray = None,
                     ce_t: np.ndarray = None, fmax: float = 100.0) -> None:
    """Time-frequency spectrogram of the 80/20 EPOP/IPOP weighted-mean LFP.

    Optionally adds an x0/CE schedule panel below when t, x0_t, and ce_t
    are provided.

    Args:
        filepaths (Any): FilePaths with figures_dir.
        params_dict (Dict): Needs TRANSIENT.
        x1 (np.ndarray): (num_cells, time) EPOP x.
        x2 (np.ndarray): (num_cells, time) IPOP x.
        fs (float): Sampling rate of x1/x2 in Hz (see synch_metrics.sample_rate).
        t (np.ndarray): Time vector matching x1/x2 (seconds). Required for param panel.
        x0_t (np.ndarray): Recorded x0 schedule over time. If None, param panel is skipped.
        ce_t (np.ndarray): Recorded CE schedule over time. If None, param panel is skipped.
        fmax (float): Upper frequency limit for display in Hz. Default 100 Hz.
    
    Returns: void
    """
    x_lfp = 0.8 * np.mean(x1, axis=0) + 0.2 * np.mean(x2, axis=0)

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

