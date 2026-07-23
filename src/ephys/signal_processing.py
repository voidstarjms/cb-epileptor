import numpy as np
from igor2 import binarywave
import matplotlib.pyplot as plt
from scipy.signal import filtfilt, butter, find_peaks, welch
import argparse
import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
import sheet_parser

sys.path.append("..")
from plotting import ephys_plots

# NOTE: The default directory must exist for this to work
DEFAULT_IN_DIR = os.path.join("..", "..", "ephys", "binaries")
DEFAULT_SHEET_PATH = os.path.join("..", "..", "ephys", "master_sheet.ods")
DEFAULT_OUT_DIR = os.path.join("..", "..", "figures", "ephys")
DRIFT_FILT_FREQ = 200
STD_THRESH = 4.75
MAX_PRE_PEP_SWEEP = -16

def analyze_single_binary(fname, pep_sweep=None, disp_sweep=-1, show=False, verbose=True, out_dir=None):
    if verbose:
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
    y_raw = np.empty((end_time - start_time, y.shape[1]))
    y_processed = np.empty((end_time - start_time, y.shape[1]))
    transients_per_sweep = np.zeros(y.shape[1], dtype=int)
    for i in range(y.shape[1]):
        y_filtered = butter_highpass_filter(y[:, i], DRIFT_FILT_FREQ, fs)
        y_raw[:, i] = y_filtered
        
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
    
    prepep_transients = 0
    postpep_transients = 0
    if disp_sweep >= 0 and verbose:
        print("Transients in specified sweep:", transients_per_sweep[disp_sweep])
    if pep_sweep != None:
        prepep_transients = np.sum(transients_per_sweep[:pep_sweep])
        postpep_transients = np.sum(transients_per_sweep[pep_sweep:])
        if verbose:
            print("Total pre-PEP transients detected:", prepep_transients)
            print("Total post-PEP transients detected:", postpep_transients)
    else:
        prepep_transients = np.sum(transients_per_sweep)
        if verbose:
            print("Total transients detected:", prepep_transients)

    # Plot the LFPs if requested
    if show == True:
        ephys_plots.plot_lfp(x, y_raw, y_processed, fname, disp_sweep, out_dir=out_dir)

    return prepep_transients, postpep_transients

def analyze_binaries(sheet_path : str, bin_dir : str, verbose=False, out_dir=None, show=False):
    sheet_df = sheet_parser.parse_sheet(sheet_path)

    subdirs = os.listdir(bin_dir)
    subdirs.sort(key=lambda x: datetime.strptime(x, '%m %d %Y'))
    i = 0
    # tuple structure is (expt_count, total_transients)
    expt_types = ["control",
                  "cbd",
                  "picro",
                  "cbd_picro",
                  "cbd_bc"]
    expt_types = ["wt_"+t for t in expt_types] +\
        ["btbr_"+t for t in expt_types]
    transients_by_expt_type = {t+"_pre": [] for t in expt_types} |\
        {t+"_post": [] for t in expt_types}

    for d in subdirs:
        subdir_path = os.path.join(bin_dir, d)
        file_list = os.listdir(subdir_path)
        file_list.sort()
        for f in file_list:
            if sheet_df["data_OK"][i] == True:
                if verbose:
                    run_id = f"{sheet_df["date"][i]} {sheet_df["run_num"][i]}"
                    print("Skipping confounded experiment " + run_id)
            else:
                fpath = os.path.join(subdir_path, f)
                pep_sweep_count = -1
                pep_sweep = "transient-1"
                while (pep_sweep_count >= MAX_PRE_PEP_SWEEP-1):
                    if sheet_df[pep_sweep][i] == "":
                        break
                    pep_sweep = "transient"+str(pep_sweep_count)
                    pep_sweep_count -= 1
                pre_counted, post_counted = analyze_single_binary(fpath, pep_sweep=pep_sweep_count, verbose=verbose)
                
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
                transients_by_expt_type[key + "_pre"].append(pre_counted)
                transients_by_expt_type[key + "_post"].append(post_counted)

                if verbose:
                    print("Expected transients:", sheet_df["total_transients"][i])
                #rel_err = abs(sheet_df["total_transients"][i] -
                #              counted) / counted
                #print(f"Relative error: {rel_err}")
                #sum_rel_err += rel_err
            i += 1
            if verbose:
                print()

    if show:
        ephys_plots.plot_scatter_columns(expt_types, transients_by_expt_type, out_dir=out_dir)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fname', type=str, default=None,
                        help='Path to ibw binary file.')
    parser.add_argument('--dirname', type=str, default=DEFAULT_IN_DIR,
                        help="""Path to directory of ibw binary files. 
Subdirectories must have naming scheme m dd yyyy.""")
    parser.add_argument('--sweep', default=0, type=int,
                        help='The 1-indexed sweep number to examine. Leave out for all sweeps')
    parser.add_argument('--sheet-path', default=DEFAULT_SHEET_PATH, type=str,
                        help="""Path to master spreadsheet file to supply ancillary data for 
binary analysis.""")
    parser.add_argument('--show', default=False, action='store_true',
                        help='Show plots after completion.')
    parser.add_argument('--out-dir', default=None, type=str,
                        help='Directory to save figures to.')
    parser.add_argument('-v', '--verbose', action='store_true',
                        default=False, help='Print binary analyses for each file')
    parser.add_argument('-q', '--quiet', action='store_true', default=False,
                        help='Print nothing. Overrides verbose')

    args = parser.parse_args()
    fname = args.fname
    dirname = args.dirname
    sheet_path = args.sheet_path

    sweep = args.sweep - 1
    show = args.show
    quiet = args.quiet
    verbose = args.verbose if not quiet else False

    if fname != None:
        analyze_single_binary(fname, disp_sweep=sweep, show=show, verbose=not quiet)
    else:
        analyze_binaries(sheet_path, dirname, show=show, verbose=verbose)
        
if __name__ == "__main__":
    main()