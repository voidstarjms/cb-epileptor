"""Synchrony measures: chi and the Kuramoto order parameter.

Provides the chi population-synchrony statistic (with Gaussian-smoothed
autocorrelation) and a spike-derived KOP that interpolates phase 2*pi*k
between consecutive spikes on a 1 ms grid.
"""
import numpy as np
import scipy
from typing import Tuple


def autocorrelate(data: np.ndarray) -> Tuple[float, np.ndarray, np.ndarray]:
    """Smooth data with a Gaussian, then return chi, autocorr, lag.

    Args:
        data (np.ndarray): (num_neurons, time) array of a state variable.

    Returns:
        Tuple[float, np.ndarray, np.ndarray]: chi (synchrony measure scalar),
        autocorr (autocorrelation of the population mean), and lag (lag values
        matching autocorr).
    """
    smoothed_data = scipy.ndimage.gaussian_filter(data, sigma=2.0)
    chi, autocorr, lag = synchrony_stats(smoothed_data)
    return chi, autocorr, lag

def synchrony_stats(data: np.ndarray, maxlags: int = 3000) -> Tuple[float, np.ndarray, np.ndarray]:
    """Compute chi (population synchrony) and the autocorrelation of the population mean.

    Args:
        data (np.ndarray): (num_neurons, time) array.
        maxlags (int): Maximal lag for autocorrelation. Defaults to 3000 ms.

    Returns:
        Tuple[float, np.ndarray, np.ndarray]: chi (sqrt of pop-mean variance
        divided by mean single-neuron variance), autocorr (autocorrelation of
        the population mean), and lag (lag values matching autocorr).
    """
    data_pop=np.mean(data, axis=0) # pop avg
    sigma_pop=np.mean(np.square(data_pop)) - np.square(np.mean(data_pop))
    sigma=np.mean(np.square(data), axis=1) - np.square(np.mean(data, axis=1))
    sigma_mean=np.mean(sigma)
    chisq=sigma_pop / sigma_mean
    chi=np.sqrt(chisq)

    mean_subtract=data_pop - np.mean(data_pop)
    autocorr=scipy.signal.correlate(mean_subtract, mean_subtract, mode='same')
    lag = scipy.signal.correlation_lags(len(mean_subtract), len(mean_subtract), mode='same')
    return chi, autocorr, lag


def KOP(neuron_idx: np.ndarray, spike_times: np.ndarray, duration: float,
        t_start: float = 0.0, time_bin_size: float = 0.001) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Kuramoto order parameter from spike-derived phases.

    Args:
        neuron_idx (np.ndarray): Per-spike neuron index.
        spike_times (np.ndarray): Per-spike times (seconds).
        duration (float): Simulation duration (seconds).
        t_start (float): Time of the first recorded sample (seconds), i.e. TRANSIENT.
        time_bin_size (float): Step of the phase time grid in seconds (default 0.001 = 1 ms).

    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray]: Z (complex mean of
        exp(i*phase) across neurons, shape (time,)), r (magnitude of Z;
        0 = no lock, 1 = full lock), and psi (angle of Z in radians, i.e.
        the collective phase).
    """
    phase = compute_phase(neuron_idx, spike_times, duration, t_start, time_bin_size)
    Z = np.mean(np.exp(1j * phase), axis=0)

    # r is the magnitude of z
    r = np.abs(Z)
    # psi is the angle of z from the real axis
    psi = np.angle(Z)
    return Z, r, psi

    
def mean_firing_rate(neuron_idx: np.ndarray, spike_times: np.ndarray,
                     num_cells: int, duration: float) -> float:
    """Mean firing rate across all neurons in Hz.

    Args:
        neuron_idx: Per-spike neuron index (unused, kept for consistent signature).
        spike_times: Per-spike times in seconds.
        num_cells: Total number of neurons in the population.
        duration: Recording duration in seconds.

    Returns:
        float: Mean firing rate in Hz (total spikes / num_cells / duration).
    """
    return len(spike_times) / (num_cells * duration)


def sample_rate(t: np.ndarray) -> float:
    """Sampling rate of recorded traces in Hz, from the recorded time vector.

    Use this instead of the sim clock: StateMonitors may record at a coarser
    dt than defaultclock.dt (see model.run_sim's low_res).

    Args:
        t: Recorded time vector in seconds.

    Returns:
        float: Sampling rate in Hz.
    """
    return float(1 / (t[1] - t[0]))


def median_frequency(x: np.ndarray, fs: float) -> float:
    """Median of per-window peak frequencies from the spectrogram.

    Computes a spectrogram of the population mean signal, finds the peak
    frequency in each time window, and returns the median
    across all windows.

    Args:
        x: (num_cells, time) array for one population.
        fs: Sampling rate of x in Hz (see sample_rate).

    Returns:
        float: Median peak frequency in Hz.
    """
    x_mean = np.mean(x, axis=0)
    f, _, Sxx = scipy.signal.spectrogram(x_mean, fs=fs, nperseg=int(fs), noverlap=int(fs) // 2)
    peak_freqs = f[1:][np.argmax(Sxx[1:], axis=0)] 
    return float(np.median(peak_freqs))


def isi_cv(neuron_idx: np.ndarray, spike_times: np.ndarray) -> float:
    """Mean coefficient of variation of inter-spike intervals across neurons.

    CV = std(ISI) / mean(ISI) per neuron, then averaged. Neurons with fewer
    than 2 spikes dont have interspike interval and are skipped. Returns NaN if no neuron has enough spikes.

    Args:
        neuron_idx: Per-spike neuron index.
        spike_times: Per-spike times in seconds.

    Returns:
        float: Mean ISI CV across neurons. 0 = perfectly regular, >1 = bursty.
    """
    cvs = []
    for idx in np.unique(neuron_idx):
        times = np.sort(spike_times[neuron_idx == idx])
        if len(times) < 2:
            continue
        intervals = np.diff(times)
        mu = np.mean(intervals)
        if mu > 0:
            cvs.append(np.std(intervals) / mu)
    return float(np.mean(cvs)) if cvs else float('nan')


def summarize_run(results: dict, num_cells: int, duration: float, t_start: float) -> dict:
    """Every summary metric for both populations of one simulation.

    Single source for the summary numbers, used by run.report_metrics and by
    sweep/build_table.py, so a metric fix here reaches both.

    Args:
        results: The 'results' dict of a saved run (t, x1, x2, spikes_n1, spikes_n2).
        num_cells: Neurons per population.
        duration: Recorded duration in seconds (SIM_DURATION, excluding the transient).
        t_start: Time of the first recorded sample in seconds (TRANSIENT).

    Returns:
        dict with chi, KOP r (full-run mean), mean firing rate, median peak
        frequency and ISI CV, each suffixed _exc (x1/spikes_n1) or _inh (x2/spikes_n2).
    """
    fs = sample_rate(results['t'])
    summary = {}
    for pop, trace, spikes in (('exc', 'x1', 'spikes_n1'), ('inh', 'x2', 'spikes_n2')):
        x = results[trace]
        idx, times = results[spikes]['i'], results[spikes]['t']
        chi, _, _ = autocorrelate(x)
        _, r, _ = KOP(idx, times, duration, t_start)
        summary[f'chi_{pop}'] = float(chi)
        summary[f'r_{pop}'] = float(np.mean(r))
        summary[f'rate_{pop}_hz'] = mean_firing_rate(idx, times, num_cells, duration)
        summary[f'median_freq_{pop}_hz'] = median_frequency(x, fs)
        summary[f'isi_cv_{pop}'] = isi_cv(idx, times)
    return summary



def compute_phase(neuron_idx: np.ndarray, spike_times: np.ndarray, duration: float,
                  t_start: float = 0.0, time_bin_size: float = 0.001) -> np.ndarray:
    """Interpolate phase on a uniform time grid: the k-th spike gets phase 2*pi*k.

    Args:
        neuron_idx (np.ndarray): Per-spike neuron index.
        spike_times (np.ndarray): Per-spike times (seconds).
        duration (float): Simulation duration (seconds).
        t_start (float): Time of the first recorded sample (seconds), i.e. TRANSIENT.
            Spike times are absolute, so the grid must start here to line up with them.
        time_bin_size (float): Step of the time grid in seconds (default 0.001 = 1 ms).

    Returns:
        np.ndarray: (num_unique_neurons, num_time_steps) phase matrix in radians.
    """
    # Uniform discritized representation of time. 
    # All oscillators will be mapped to this scale
    # these are also the time steps of theta in the KOP
    time_grid = np.arange(t_start, t_start+duration+time_bin_size, time_bin_size)

    # find unique neuron indices to iterate through
    unique_idx = np.unique(neuron_idx)

    # matrix 'theta' of the KOP
    # dim: num_oscillators x num_time_steps
    phase_matrix = []

    for idx in unique_idx:
        # find time of spike for each neuron oscillator
        mask = np.where(neuron_idx == idx)
        x_coord = spike_times[mask]
        y_coord =  2 * np.pi *  np.arange(0, len(x_coord))

        # 'map' spikes to multiples of unit circle by setting 
        # spike points as 2pi*k and 2pi*(k+1) 
        # and interpolating time steps between them. 
        theta = np.interp(time_grid, x_coord, y_coord)
        phase_matrix.append(theta)

    return np.array(phase_matrix)
