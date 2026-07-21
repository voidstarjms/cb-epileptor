import numpy as np
from igor2 import binarywave
import matplotlib.pyplot as plt
from scipy.signal import filtfilt, butter, find_peaks
from sklearn.decomposition import FastICA
from sklearn.preprocessing import StandardScaler
import argparse
import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from master_sheet_processing import sheet_parser

DRIFT_FILT_FREQ = 200
STD_THRESH = 4.75

def analyze_single_binary(fname, sweep=-1, show=False):
    print(f"Loading binary {Path(fname).stem}")
    wave = binarywave.load(fname)
    w = wave['wave']

    # Data
    y = w['wData'] # 1D numpy array of your waveform values

    fs = 10000

    def butter_highpass_filter(data, cutoff, fs, order=2):
        nyq = 0.5 * fs
        high = cutoff / nyq
        b, a = butter(order, high, btype='high')
        return filtfilt(b, a, data)
    def butter_lowpass_filter(data, cutoff, fs, order=2):
        nyq = 0.5 * fs
        low = cutoff / nyq
        b, a = butter(order, low, btype='low')
        return filtfilt(b, a, data)
    
    start_time = 0
    end_time = 200000
    y_processed = np.empty((end_time - start_time, y.shape[1]))
    transients_per_sweep = np.zeros(y.shape[1], dtype=int)
    for i in range(y.shape[1]):
        y_filtered = butter_highpass_filter(y[:, i], DRIFT_FILT_FREQ, fs)
        mean = np.mean(y_filtered)
        stdev = np.std(y_filtered)
        y_filtered = y_filtered[start_time:end_time] - mean
        y_thresholded = np.where(y_filtered > STD_THRESH * stdev, y_filtered, 0)
        transients_per_sweep[i] = find_peaks(y_thresholded)[0].shape[0]
        y_processed[:, i] = y_thresholded

    # Reconstruct x axis
    header = w['wave_header']
    delta  = header['sfA'][0]   # spacing between points
    offset = header['sfB'][0]   # x value of first point

    x = offset + delta * np.arange(len(y)) + start_time * delta
    # ica = FastICA(n_components=2)
    # ica.fit(y.T)
    # print(ica.mixing_.shape)
    # S_estimated = ica.transform(y.T[0].reshape(1, -1))
    
    if sweep >= 0:
        print("Transients in this sweep:", transients_per_sweep[sweep])
    total_transients = np.sum(transients_per_sweep)
    print("Total transients detected:", total_transients)

    if show == True:
        # Plot all LFPs overlaid
        if sweep == -1:
            fig = plt.figure()
            ax = fig.add_axes((0.1, 0.1, 0.8, 0.8))
            fig.set_size_inches(7, 6)
            fname_parts = fname.split(os.sep)[-1].split(".")[0].split("_")
            sweep_num = f"Sweep {sweep+1}" if sweep >= 0 else "All Sweeps"
            fdate = f"{fname_parts[-3]}-{fname_parts[-2]}-{fname_parts[-1]} {fname_parts[0]} {sweep_num}"
            fig.suptitle(f"{fdate}")
            for i in range(y.shape[1]):
                ax.plot(x, butter_highpass_filter(y[:, i], DRIFT_FILT_FREQ, fs), alpha=0.4)

            ax.set_xlabel("Time (s)")
            ax.set_ylabel("Raw LFP (V)")
        # PLot a single LFP and the transients detected from it
        else:
            fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
            fig.set_size_inches(14, 6)
            fname_parts = fname.split(os.sep)[-1].split(".")[0].split("_")
            sweep_num = f"Sweep {sweep+1}" if sweep >= 0 else "All Sweeps"
            fdate = f"{fname_parts[-3]}-{fname_parts[-2]}-{fname_parts[-1]} {fname_parts[0]} {sweep_num}"
            fig.suptitle(f"{fdate}")
            ax1.plot(x, butter_highpass_filter(y[:, sweep], DRIFT_FILT_FREQ, fs))
            ax2.plot(x, y_processed[:, sweep])

            ax2.set_xlabel("Time (s)")
            ax1.set_ylabel("Raw LFP (V)")
            ax2.set_ylabel("Processed LFP\n(V relative to mean)")

        #plt.xlabel(w['dimension_units'].decode().strip('\x00'))
        #plt.ylabel(w['data_units'].decode().strip('\x00'))

        plt.show()

    return total_transients

def compare_sheet_to_binaries(sheet_path : str, bin_dir : str):
    sheet_df = sheet_parser.parse_sheet(sheet_path)

    subdirs = os.listdir(bin_dir)
    subdirs.sort(key=lambda x: datetime.strptime(x, '%m %d %Y'))
    i = 0
    sum_rel_err = 0
    # tuple structure is (expt_count, total_transients)
    transients_by_expt_type = {"wt_control" : [],
                    "wt_cbd" : [],
                    "wt_picro" : [],
                    "wt_cbd_picro" : [],
                    "wt_cbd_bc" : [],
                    "btbr_control" : [],
                    "btbr_cbd" : [],
                    "btbr_picro" : [],
                    "btbr_cbd_picro" : [],
                    "btbr_cbd_bc" : []}
    for d in subdirs:
        subdir_path = os.path.join(bin_dir, d)
        file_list = os.listdir(subdir_path)
        file_list.sort()
        for f in file_list:
            if sheet_df["data_OK"][i] == True:
                run_id = f"{sheet_df["date"][i]} {sheet_df["run_num"][i]}"
                print("Skipping confounded experiment " + run_id)
            else:
                fpath = os.path.join(subdir_path, f)
                counted_transients = analyze_single_binary(fpath)
                
                # Add transient counts to appropriate tuple
                has_cbd = sheet_df["has_cbd"][i]
                has_picro = sheet_df["has_picro"][i]
                has_bc = sheet_df["has_bc"][i]
                key = ""
                key += "cbd" if has_cbd else ""
                key += (("_" if key != "" else "") + "picro") if has_picro else ""
                key += "_bc" if has_bc else ""
                key = key if key != "" else "control"
                key = ("btbr_" if sheet_df["genotype"][i] == "BTBR" else "wt_") + key
                transients_by_expt_type[key].append(counted_transients)
                
                print("Expected transients:", sheet_df["total_transients"][i])
                rel_err = abs(sheet_df["total_transients"][i] -
                              counted_transients) / counted_transients
                print(f"Relative error: {rel_err}")
                sum_rel_err += rel_err
            i += 1
            print()
    
    print("-"*40)
    print(f"Mean relative error: {sum_rel_err / i}\n")
    print("Mean and CV of transient counts by experiment type\n")
    coeffs_of_var = np.zeros(len(transients_by_expt_type))
    label_width = max(len(k) for k in transients_by_expt_type.keys())
    # Print header
    print(f"{"Run type":<{label_width}}{"Mean":>10}{"CV":>10}{"N":>10}")
    # Compute and print mean and coefficient of variation for each expt type
    for i, (k, v) in enumerate(transients_by_expt_type.items()):
        run_count = len(v)
        transient_arr = np.asarray(v)
        mean = 0
        std = 0
        cv = 0
        # Avoid divide by zero if no runs of given type
        if run_count != 0:
            mean = np.mean(transient_arr)
            std = np.std(transient_arr)
            cv = std / mean if mean != 0 else np.nan
        coeffs_of_var[i] = cv
        print(f"{k:<{label_width}}{mean:>10.3g}{cv:>10.3g}{run_count:>10}")
    print(f"\nMean coefficient of variation: {np.mean(coeffs_of_var)}")



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fname', type=str, default=None)
    parser.add_argument('--dirname', type=str, default=None)
    parser.add_argument('--sweep', default=0, type=int)
    parser.add_argument('--sheet-path', default=None,
                        type=str)
    parser.add_argument('--show', default=False, action='store_true')

    args = parser.parse_args()
    fname = args.fname
    dirname = args.dirname
    sheet_path = args.sheet_path

    # Argument sanity checking
    if fname == None and dirname == None:
        print("You must specify either a file or directory for input")
        sys.exit(1)
    if fname != None and dirname != None:
        print("You cannot specify both a file and a directory")
        sys.exit(2)
    if dirname != None and sheet_path == None:
        print("You must provide a spreadsheet to compare against a directory's contents. Please pass --sheet-path [path]")
        sys.exit(3)    
    
    sweep = args.sweep - 1
    show = args.show

    if fname != None:
        analyze_single_binary(fname, sweep=sweep, show=show)
    else:
        compare_sheet_to_binaries(sheet_path, dirname)
        
if __name__ == "__main__":
    main()