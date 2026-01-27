import numpy as np
import scipy

def autocorelate(data):
    # apply gaussian smoothing
    smoothed_data = scipy.ndimage.gaussian_filter(data, sigma=2.0)
    chi, autocorr = synchrony_stats(smoothed_data)
    return chi, autocorr

def synchrony_stats(data, maxlags=3000):
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
    autocorr=scipy.signal.correlate(mean_subtract, mean_subtract, mode='valid')
    return chi, autocorr


def KOP():

    pass