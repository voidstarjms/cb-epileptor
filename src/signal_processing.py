import numpy as np
from igor2 import binarywave
import matplotlib.pyplot as plt
from scipy.signal import filtfilt, butter, find_peaks
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--fname')
parser.add_argument('--cutoff', type=float)

args = parser.parse_args()

wave = binarywave.load(args.fname)
w = wave['wave']

# Data
y = w['wData'] # 1D numpy array of your waveform values

fs = 10000

def butter_highpass_filter(data, cutoff, fs, order=2):
    nyq = 0.5 * fs
    high = cutoff / nyq
    b, a = butter(order, high, btype='high')
    return filtfilt(b, a, data)

cutoff = args.cutoff

end_time = 180000
y_processed = np.empty((end_time, y.shape[1]))
transient_total = 0
for i in range(y.shape[1]):
    y_filtered = butter_highpass_filter(y[:end_time, i], cutoff, fs)
    y_rectified = np.where(y_filtered < 0, 0, y_filtered)
    peak_mask = np.full(end_time, False)
    peak_mask[find_peaks(y_rectified)[0]] = True
    mean_peak_height = np.mean(y_rectified[:end_time], where=peak_mask)
    print(mean_peak_height)
    y_thresholded = np.where(np.abs(y_rectified) < 2 * mean_peak_height, 0, y_rectified)
    transient_total += find_peaks(y_thresholded)[0].shape[0]
    y_processed[:, i] = y_thresholded
y = y_processed
print(transient_total)

# Reconstruct x axis
header = w['wave_header']
delta  = header['sfA'][0]   # spacing between points
offset = header['sfB'][0]   # x value of first point

x = offset + delta * np.arange(len(y))

plt.plot(x, y)
plt.xlabel(w['dimension_units'].decode().strip('\x00'))
plt.ylabel(w['data_units'].decode().strip('\x00'))
plt.show()
