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
MAX_PRE_PEP_SWEEP = sheet_parser.MAX_PRE_PEP_SWEEP
MAX_POST_PEP_SWEEP = sheet_parser.MAX_POST_PEP_SWEEP

def analyze_single_binary(fname, df=None, entry_idx=None, df_start_pos=None, disp_sweep=-1, show=False, verbose=True, out_dir=None):
    """Load and count the transients in a single ibw binary file.

    Args:
        fname (str): The path of the input ibw file to load.
        df (DataFrame): The DataFrame containing metadata from the master sheet
        entry_idx (int): The index of the experiment's entry in the DataFrame. Must be set if df is not None.
        df_start_pos (int): The number of sweeps before PEP introduction. Must be set if df is not None.
        disp_sweep (int): The zero-indexed sweep from the input file to plot. -1 plots all sweeps together.
        show (bool): Show the plot.
        verbose (bool): Print extra output info.
        out_dir (str): The directory to save figures to. Saves to DEFAULT_OUT_DIR if None.
    
    Returns:
        tuple (prepep_transients, postpep_transients): total number of transients counted before and after PEP
            introduction. prepep_transients is all transients and postpep_transients is 0 if df is not supplied.
    """
    if type(df) is not type(None) and (df_start_pos == None or entry_idx == None):
        print("analyze_single_binary: You must specify a row and transient starting column to analyze using a spreadsheet")
        sys.exit(1)
    scan_for_nan = type(df) is not type(None)
    # Compute the sweep that PEP is introduced
    pep_sweep = -df_start_pos if df_start_pos != None else None

    # Load the binary
    if verbose:
        print(f"Loading binary {Path(fname).stem}")
    wave = binarywave.load(fname)
    w = wave['wave']
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
        # Check for NaN sweep, continue if so
        if scan_for_nan:
            df_sweep_name = "transient" + ("+" if df_start_pos >= 0 else "") + str(df_start_pos)
            df_start_pos += 1
            if df[df_sweep_name][entry_idx] == np.nan:
                 continue

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
    if show:
        ephys_plots.plot_lfp(x, y_raw, y_processed, fname, disp_sweep, out_dir=out_dir)

    return prepep_transients, postpep_transients

def _find_first_transient_cell(sheet_df : pd.DataFrame, idx : int):
    # Find the starting sweep in the DataFrame
    start_sweep_count = -1
    while (start_sweep_count >= MAX_PRE_PEP_SWEEP-1):
        start_sweep = "transient"+str(start_sweep_count)
        if sheet_df[start_sweep][idx] == "":
            start_sweep_count += 1
            break
        start_sweep_count -= 1

    return start_sweep_count

def _derive_experiment_key(sheet_df : pd.DataFrame, idx : int):
    has_cbd = sheet_df["has_cbd"][idx]
    has_picro = sheet_df["has_picro"][idx]
    has_bc = sheet_df["has_bc"][idx]

    key = ""
    key += "cbd" if has_cbd else ""
    key += (("_" if key != "" else "") + "picro") if has_picro else ""
    key += "_bc" if has_bc else ""
    key = key if key != "" else "control"
    key = ("btbr_" if sheet_df["genotype"][idx] == "BTBR" else "wt_") + key

    return key

def _get_manual_transient_count_single(df : pd.DataFrame, idx : int):
    expt_type = _derive_experiment_key(df, idx)
    sweep_idx = _find_first_transient_cell(df, idx)
    
    sweep_name = "transient"+("" if sweep_idx < 0 else "+")+str(sweep_idx)
    prepep_transients = 0
    postpep_transients = 0
    while (sweep_idx < MAX_POST_PEP_SWEEP and df[sweep_name][idx] != ""):
        if not np.isnan(df[sweep_name][idx]):
            if sweep_idx < 0:
                prepep_transients += df[sweep_name][idx]
            else:
                postpep_transients += df[sweep_name][idx]
        sweep_idx += 1
        sweep_name = "transient"+("" if sweep_idx < 0 else "+")+str(sweep_idx)

    return expt_type, prepep_transients, postpep_transients

def _init_transient_dict():
    # Create set of lists to hold transients counts for each experiment type
    expt_types = ["control",
                  "cbd",
                  "picro",
                  "cbd_picro",
                  "cbd_bc"]
    expt_types = ["wt_"+t for t in expt_types] +\
        ["btbr_"+t for t in expt_types]
    transients_by_expt_type = {t+"_pre": [] for t in expt_types} |\
        {t+"_post": [] for t in expt_types}
    
    return expt_types, transients_by_expt_type

def get_manual_transient_count(df : pd.DataFrame):
    row_count = len(df)
    expt_types, transients = _init_transient_dict()
    for i in range(row_count):
        if df["data_OK"][i] == True:
            continue
        expt, prepep, postpep = _get_manual_transient_count_single(df, i)
        transients[expt + "_pre"].append(prepep)
        transients[expt + "_post"].append(postpep)

    return expt_types, transients

def analyze_binaries(sheet_df : pd.DataFrame, bin_dir : str, verbose=False):
    """Analyze a directory of binary files using a master spreadsheet as metadata.
        Binaries must be organized into subdirectories by date of experiment, and
        subdirectories must be named in the format '[m]m dd yyyy'. Binaries must be
        named starting with 'R[run number]_'.

        Args:
            sheet_df (DataFrame): The processed DataFrame derived from the
                master spreadsheet.
            bin_dir (str): The path to the topmost binary directory.
            verbose (bool): Print extra output info.
        
        Returns:
            tuple (expt_types, transients_by_expt_type):
                expt_types: list of string keys for each experiment types.
                transients_by_expt_type: dict of lists of transient counts
                    in each experiment for each type.
    """
    # Load all subdirectories and sort by date
    subdirs = os.listdir(bin_dir)
    subdirs.sort(key=lambda x: datetime.strptime(x, '%m %d %Y'))

    expt_types, transients_by_expt_type = _init_transient_dict()
    # Iterate over all subdirectories
    i = 0 # Entry index in the DataFrame
    for d in subdirs:
        subdir_path = os.path.join(bin_dir, d)
        file_list = os.listdir(subdir_path)
        file_list.sort()

        # Iterate over all files in the subdirectory
        for f in file_list:
            start_sweep_count = _find_first_transient_cell(sheet_df, i)
            
            if sheet_df["data_OK"][i] == True:
                if verbose:
                    run_id = f"{sheet_df["date"][i]} {sheet_df["run_num"][i]}"
                    print("Skipping confounded experiment " + run_id)
            else:
                fpath = os.path.join(subdir_path, f)
                pre_counted, post_counted = analyze_single_binary(fpath, df=sheet_df, entry_idx=i,
                                                                    df_start_pos=start_sweep_count-1,
                                                                    verbose=verbose, show=False)

                # Add transient counts to appropriate lists
                key = _derive_experiment_key(sheet_df, i)
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

    return expt_types, transients_by_expt_type

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--fname', type=str, default=None,
                        help='Path to ibw binary file.')
    parser.add_argument('--dirname', type=str, default=DEFAULT_IN_DIR,
                        help="""Path to directory of ibw binary files. 
Subdirectories must have naming scheme m dd yyyy.""")
    parser.add_argument('--sweep', default=0, type=int,
                        help='The 1-indexed sweep number to examine. Leave out for all sweeps')
    parser.add_argument('--sheet-path', default=DEFAULT_SHEET_PATH, type=str,
                        help="""Path to master spreadsheet file to supply ancillary data for binary analysis.""")
    parser.add_argument('--show', default=False, action='store_true',
                        help='Show plots after completion.')
    parser.add_argument('--out-dir', default=DEFAULT_OUT_DIR, type=str,
                        help='Directory to save figures to.')
    parser.add_argument('-v', '--verbose', action='store_true',
                        default=False, help='Print binary analyses for each file.')
    parser.add_argument('-q', '--quiet', action='store_true', default=False,
                        help='Print nothing. Overrides verbose.')
    parser.add_argument('--mode', type=str, default='scatter_col_auto',
                        help="""Mode for multi-binary analysis plots.
    \nscatter_col_auto: Plot scatter column plot of pre- and post-PEP transient counts by
    experiment type, using std-threshold signal analysis.
    \nscatter_col_man: Plot scatter column plot of pre- and post-PEP transient counts by
    experiment type, using manual transient tallies from master sheet.
    \nmean_cv: Plot two-row bar graph of mean and CV of transient counts by experiment
    type, plotting both pre- and post-PEP.
    \nauto_v_man: Plot a scatter plot of the automatic transient tally against the
    manual one.""")

    args = parser.parse_args()
    fname = args.fname
    dirname = args.dirname
    sheet_path = args.sheet_path
    out_dir = args.out_dir
    sweep = args.sweep - 1
    show = args.show
    quiet = args.quiet
    verbose = args.verbose if not quiet else False
    mode = args.mode

    if fname != None:
        analyze_single_binary(fname, disp_sweep=sweep, show=show, verbose=not quiet, out_dir=out_dir)
    else:
        sheet_df = sheet_parser.parse_sheet(sheet_path)
        if "auto" in mode:
            expt_types, transients_auto = analyze_binaries(sheet_df, dirname, verbose=verbose)
        if "man" in mode:
            expt_types, transients_man = get_manual_transient_count(sheet_df)
        match mode:
            case 'scatter_col_auto':
                ephys_plots.plot_scatter_columns(expt_types, transients_auto,
                                                 show=show, out_dir=DEFAULT_OUT_DIR, out_suffix="_auto",
                                                 title="Transients by Experiment (Automated Counting)")
            case 'scatter_col_man':
                ephys_plots.plot_scatter_columns(expt_types, transients_man,
                                                 show=show, out_dir=DEFAULT_OUT_DIR, out_suffix="_man",
                                                 title="Transients by Experiment (Manual Counting)")
            case 'mean_cv':
                pass
            case 'auto_v_man':
                pass
        
if __name__ == "__main__":
    main()