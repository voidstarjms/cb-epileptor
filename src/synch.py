
import numpy as np
import scipy
from typing import Tuple


def autocorrelate(data: np.ndarray) -> Tuple[float, np.ndarray, np.ndarray]:
    # apply gaussian smoothing
    smoothed_data = scipy.ndimage.gaussian_filter(data, sigma=2.0)
    chi, autocorr, lag = synchrony_stats(smoothed_data)
    return chi, autocorr, lag

def synchrony_stats(data: np.ndarray, maxlags: int = 3000) -> Tuple[float, np.ndarray, np.ndarray]:
    '''
    Synchrony measures

    Parameters
    ==========
      data: numneuron x time
      dt: time spacing
      maxlags: maximal lag for autocorrelation, default=3000 ms

    Returns
    =======
      chi: synchrony measure
      autocorr: autocorrelation of population avg \bar{data}(t)
    '''
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
    phase = compute_phase(neuron_idx, spike_times, duration)
    Z = np.mean(np.exp(1j * phase), axis=0)

    # r is the magnitude of z
    r = np.abs(Z)
    # psi is the angle of z from the real axis
    psi = np.angle(Z)
    return Z, r, psi

    
def compute_phase(neuron_idx: np.ndarray, spike_times: np.ndarray, duration: float) -> np.ndarray:
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
