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


def KOP(neuron_idx: np.ndarray, spike_times: np.ndarray, duration: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Kuramoto order parameter from spike-derived phases.

    Args:
        neuron_idx (np.ndarray): Per-spike neuron index.
        spike_times (np.ndarray): Per-spike times (seconds).
        duration (float): Simulation duration (seconds).

    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray]: Z (complex mean of
        exp(i*phase) across neurons, shape (time,)), r (magnitude of Z;
        0 = no lock, 1 = full lock), and psi (angle of Z in radians, i.e.
        the collective phase).
    """
    phase = compute_phase(neuron_idx, spike_times, duration)
    Z = np.mean(np.exp(1j * phase), axis=0)

    # r is the magnitude of z
    r = np.abs(Z)
    # psi is the angle of z from the real axis
    psi = np.angle(Z)
    return Z, r, psi

    
def _norm_firing_rate(neuron_idx: np.ndarray, spike_times: np.ndarray, duration: float) -> Tuple[float, float]:
    """Min-max normalised mean firing rate.

    Min is anchored at 0 (silence); max is the observed peak single-neuron rate.

    Returns:
        Tuple[float, float]: (norm_FR in [0, 1], mean_FR in Hz).
    """
    unique, counts = np.unique(neuron_idx, return_counts=True)
    n_neurons = len(unique)
    mean_FR = len(spike_times) / (n_neurons * duration)
    max_FR = counts.max() / duration
    norm_FR = mean_FR / max_FR if max_FR > 0 else 0.0
    return norm_FR, mean_FR


def seizure_kop(
    neuron_idx: np.ndarray,
    spike_times: np.ndarray,
    duration: float,
) -> Tuple[float, float, float]:
    """Seizure score: rate-weighted Kuramoto order parameter.

    Multiplies mean KOP magnitude by min-max normalised firing rate so that
    periods of dense continuous firing (seizure) score high and sparse-but-
    synchronised periods (post-ictal) score low.

    Args:
        neuron_idx: Per-spike neuron index.
        spike_times: Per-spike times (seconds).
        duration: Recording duration in seconds (after transient).

    Returns:
        Tuple[float, float, float]: (score in [0, 1], mean KOP r, norm_FR).
    """
    _, r_ts, _ = KOP(neuron_idx, spike_times, duration)
    r = float(np.mean(r_ts))
    norm_FR, _ = _norm_firing_rate(neuron_idx, spike_times, duration)
    score = r * norm_FR
    return score, r, norm_FR


def seizure_kop_combined(
    neuron_idx: np.ndarray,
    spike_times: np.ndarray,
    duration: float,
    alpha: float = 0.5,
    beta: float = 0.5,
) -> Tuple[float, float, float]:
    """Seizure score: linear combination of KOP and normalised firing rate.

    Score is in [0, 1] when alpha + beta == 1. Increasing alpha weights
    synchrony more; increasing beta weights firing density more.

    Args:
        neuron_idx: Per-spike neuron index.
        spike_times: Per-spike times (seconds).
        duration: Recording duration in seconds (after transient).
        alpha: Weight for mean KOP r. Defaults to 0.5.
        beta: Weight for normalised firing rate. Defaults to 0.5.

    Returns:
        Tuple[float, float, float]: (score, mean KOP r, norm_FR).
    """
    _, r_ts, _ = KOP(neuron_idx, spike_times, duration)
    r = float(np.mean(r_ts))
    norm_FR, _ = _norm_firing_rate(neuron_idx, spike_times, duration)
    score = alpha * r + beta * norm_FR
    return score, r, norm_FR


def compute_phase(neuron_idx: np.ndarray, spike_times: np.ndarray, duration: float) -> np.ndarray:
    """Interpolate phase on a 1 ms grid: the k-th spike gets phase 2*pi*k.

    Args:
        neuron_idx (np.ndarray): Per-spike neuron index.
        spike_times (np.ndarray): Per-spike times (seconds).
        duration (float): Simulation duration (seconds).

    Returns:
        np.ndarray: (num_unique_neurons, num_time_steps) phase matrix in radians.
    """
    time_bin_size = 0.001 # 0.001 of a second = millisecond
    # Uniform discritized representation of time. 
    # All oscillators will be mapped to this scale
    # these are also the time steps of theta in the KOP
    time_grid = np.arange(0, duration+time_bin_size, time_bin_size)

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
