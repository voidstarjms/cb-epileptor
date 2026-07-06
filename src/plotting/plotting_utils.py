import numpy as np

def bin_data(data : np.ndarray, bin_size : int):
    """Bin a 1D array.

    Args:
        data (ndarray): a 1D numpy array
        bin_size (int): the size of the bins
    """
    bin_means = np.zeros(int(len(data) / bin_size))
    for i in range(len(bin_means)):
        bin_means[i] = data[i*bin_size:(i+1)*bin_size].mean()
    for i in range(0, len(data)):
        data[i] = bin_means[int(i / bin_size)]
    return data