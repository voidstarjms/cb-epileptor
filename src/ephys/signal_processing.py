import numpy as np
from igor2 import binarywave
import matplotlib.pyplot as plt
from scipy.stats import median_abs_deviation, linregress
from scipy.signal import filtfilt, butter, find_peaks
import argparse
import os
import sys
import pandas as pd
from pathlib import Path
from datetime import datetime
import sheet_parser
import ephys_util as eh

sys.path.append("..")
from plotting import ephys_plots
from plotting import analysis_plots

# NOTE: The default directory must exist for this to work
DEFAULT_IN_DIR = os.path.join("..", "..", "ephys", "binaries")
DEFAULT_SHEET_PATH = os.path.join("..", "..", "ephys", "master_sheet.ods")
DEFAULT_OUT_DIR = os.path.join("..", "..", "figures", "ephys")
DEFAULT_EXPERIMENT_FILE = "expt_types.csv"
DRIFT_FILT_FREQ = 120
STD_THRESH = 7.6
EPHYS_FS = 10000
MAX_PRE_PEP_SWEEP = sheet_parser.MAX_PRE_PEP_SWEEP
MAX_POST_PEP_SWEEP = sheet_parser.MAX_POST_PEP_SWEEP

def butter_highpass_filter(data, cutoff, fs, order=2):
    nyq = 0.5 * fs
    high = cutoff / nyq
    b, a = butter(order, high, btype='high')
    return filtfilt(b, a, data)
def butter_bandstop_filter(data, lowcut, highcut, fs, order=2):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='bandstop')
    return filtfilt(b, a, data)

def get_lfp_list(fname, df=None, entry_idx=None, df_start_pos=None):
    """"""
    if type(df) is not type(None) and (entry_idx == None or df_start_pos == None):
        print("get_lfp_list: You must specify a row and transient starting column to analyze using a spreadsheet")
        sys.exit(1)
    scan_for_nan = type(df) is not type(None)

    wave = binarywave.load(fname)
    w = wave['wave']
    y = w['wData'] # numpy array of waveform values

    lfp_list = []
    for i in range(y.shape[1]):
        if scan_for_nan:
            df_sweep_name = "transient" + ("+" if df_start_pos >= 0 else "") + str(df_start_pos)
            df_start_pos += 1
            if np.isnan(float(df[df_sweep_name][entry_idx])):
                lfp_list.append(np.nan)
                continue
        lfp_list.append(butter_highpass_filter(y[:, i], 0.25, EPHYS_FS))

    return lfp_list

def _analyze_single_binary(fname, df=None, entry_idx=None, df_start_pos=None,
                            detect_freq = DRIFT_FILT_FREQ, detect_thresh = STD_THRESH,
                            verbose=True):
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

    # Load the binary
    if verbose:
        print(f"Loading binary {Path(fname).stem}")
    wave = binarywave.load(fname)
    w = wave['wave']
    y = w['wData'] # numpy array of waveform values
        
    start_time = 0
    end_time = 200000
    y_raw = np.empty((end_time - start_time, y.shape[1]))
    y_filtered = np.empty((end_time - start_time, y.shape[1]))
    y_processed = np.empty((end_time - start_time, y.shape[1]))
    thresh = np.empty(y.shape[1])
    transients_per_sweep = []
    if verbose:
        print("Sweeps in binary:", y.shape[1])
    for i in range(y.shape[1]):
        # Check for NaN sweep, continue if so
        if scan_for_nan:
            df_sweep_name = "transient" + ("+" if df_start_pos >= 0 else "") + str(df_start_pos)
            df_start_pos += 1
            if np.isnan(float(df[df_sweep_name][entry_idx])):
                transients_per_sweep.append(np.nan)
                continue

        y_filtered[:, i] = butter_highpass_filter(y[start_time:end_time, i], detect_freq, EPHYS_FS)
        y_raw[:, i] = y[:, i]
        
        median = np.median(y_filtered[:, i])
        mad = median_abs_deviation(y_filtered[:, i])

        normalized_y = y_filtered[:, i] - median
        thresh[i] = -detect_thresh * mad
        y_thresholded = np.where(normalized_y < thresh[i], normalized_y, 0)
        transients_per_sweep.append(find_peaks(-y_thresholded)[0].shape[0])
        y_processed[:, i] = y_thresholded

    # Reconstruct x axis
    header = w['wave_header']
    delta  = header['sfA'][0]   # spacing between points
    offset = header['sfB'][0]   # x value of first point

    x = offset + delta * np.arange(len(y)) + start_time * delta

    return x, y_raw, y_filtered, y_processed, transients_per_sweep, thresh

def get_binary_transients(fname, df=None, entry_idx=None, df_start_pos=None, disp_sweep=-1,
                          detect_freq = DRIFT_FILT_FREQ, detect_thresh = STD_THRESH, show=False,
                                   verbose=True, out_dir=None, scan_for_nan=True):
    """Get transient counts from a provided binary file. Supports collation with a spreadsheet for
        analysis with metadata and binning of transient counts (e.g. for pre-post analysis of
        chemical introduction).

    Args:
        fname (str): Path to binary file.
        df (DataFrame): The DataFrame containing metadata from the master sheet.
        entry_idx (int): The index of the experiment's entry in the DataFrame. Must be set if df is not None.
        df_start_pos (int): The index of the experiment's entry in the DataFrame. Must be set if df is not None.
        disp_sweep (int): The zero-indexed sweep from the input file to plot. -1 plots all sweeps together.
        show (bool): Show the plot.
        verbose (bool): Print extra output info.
        out_dir (str): The directory to save figures to. Saves to DEFAULT_OUT_DIR if None.

    Returns:
        tuple[np.array]:
            prepep_transients: Array of transient counts for pre-PEP sweeps if split point
                is specified, otherwise array of transient counts for all sweeps.
            postpep_transients: Array of transient counts for post-PEP sweeps. Unused if
                split point is not specified.
    """
    x, y_raw, y_filtered, y_processed, transients_per_sweep, std_thresh = _analyze_single_binary(fname, df=df, entry_idx=entry_idx,
                                                                         df_start_pos=df_start_pos,
                                                                         detect_freq=detect_freq, detect_thresh=detect_thresh,
                                                                         verbose=verbose)

    # Compute the sweep that PEP is introduced
    pep_sweep = -df_start_pos if df_start_pos != None else None

    # Remove nan entries and decrement PEP start position as necessary
    if pep_sweep != None and scan_for_nan == True:
        transients_nonan = []
        i = 0
        for ent in transients_per_sweep:
            if np.isnan(ent):
                if pep_sweep > i:
                    pep_sweep -= 1
            else:
                transients_nonan.append(ent)
                i += 1
        transients_per_sweep = transients_nonan

    prepep_transients = 0
    postpep_transients = 0
    if disp_sweep >= 0 and verbose:
        print("Transients in specified sweep:", transients_per_sweep[disp_sweep])
    if pep_sweep != None:
        prepep_transients = transients_per_sweep[:pep_sweep]
        postpep_transients = transients_per_sweep[pep_sweep:]
    else:
        prepep_transients = transients_per_sweep
        if verbose:
            print("Total transients detected:", prepep_transients)

    # Plot the LFPs if requested
    if show:
        ephys_plots.plot_lfp(x, y_raw, y_filtered, y_processed,
                             std_thresh, fname, disp_sweep, out_dir=out_dir)

    return prepep_transients, postpep_transients

def _find_first_transient_cell_core(df_ent : pd.Series):
    # Find the starting sweep in the DataFrame
    start_sweep_count = -1
    while (start_sweep_count >= MAX_PRE_PEP_SWEEP-1):
        start_sweep = "transient"+str(start_sweep_count)
        if df_ent[start_sweep] == "":
            # Increment back up to the last cell with content
            start_sweep_count += 1
            break
        start_sweep_count -= 1

    return start_sweep_count

def _find_first_transient_cell(sheet_df : pd.DataFrame, idx : int):
    """Compute the first used transient cell.
    
    Args:
        sheet_df (DataFrame): The DataFrame containing metadata from the master sheet.
        idx (int): Index of the entry in sheet_df.
    
    Returns:
        start_sweep_count (int): how many sweeps before PEP+0 the first sweep occurs at.
    """
    return _find_first_transient_cell_core(sheet_df.iloc[idx])

def _find_first_transient_cell_by_file(sheet_df : pd.DataFrame, fname : str):
    """"""
    file_basename = os.path.splitext(os.path.basename(fname))[0]
    fname_parts = file_basename.split("_")
    fdate = datetime(int(fname_parts[-1]), int(fname_parts[-3]),
                        int(fname_parts[-2])).date()
    run_num = int(fname_parts[0][1:])
    df_slice = sheet_df[(sheet_df["date"] == fdate) & (sheet_df["run_num"] == run_num)]
    return (_find_first_transient_cell_core(df_slice.squeeze(axis=0)), df_slice.index[0])

def _derive_experiment_key(sheet_df : pd.DataFrame, idx : int):
    """Derive the key string for an experiment type.
    
    Args:
        sheet_df (DataFrame): The DataFrame containing metadata from the master sheet.
        idx (int): Index of the entry in sheet_df.
    
    Returns:
        key (str): The key string for an experiment
    """

    has_cbd = sheet_df["has_cbd"][idx]
    has_picro = sheet_df["has_picro"][idx]
    has_bc = sheet_df["has_bc"][idx]

    key = ""
    if has_cbd:
        key += "cbd_"
        #Distinguish CBD concentrations by date
        if sheet_df["date"][idx] > datetime(2026, 3, 5).date():
            key += "10"
        else:
            key += "100"
        key += "uM"
    key += (("_" if key != "" else "") + "picro") if has_picro else ""
    key += "_bc" if has_bc else ""
    key = key if key != "" else "control"
    key = ("btbr_" if sheet_df["genotype"][idx] == "BTBR" else "wt_") + key

    return key

def get_sweeps_by_expt_type(sheet_df : pd.DataFrame, bin_dir : str, type_key : str, verbose : bool = False):
    """Get the transient counts for all sweeps in all runs of a given type.

    Args:
        sheet_df (DataFrame): DataFrame containing metadata from the master sheet.
        bin_dir (str): Path to the directory containing experiment binaries.
        type_key (str): Key string for an experiment type
        verbose (bool): Print extra output info
    
    Returns:
        tuple[np.array]:
            prepep_array (np.array): Array of transient counts for each pre-PEP sweep
                in each run of the given type, in order.
            postpep_array (np.array): Array of transient counts for each post-PEP sweep
                in each run of the given type, in order.
    """
    subdirs = os.listdir(bin_dir)
    subdirs.sort(key=lambda x: datetime.strptime(x, '%m %d %Y'))

    # Iterate over all subdirectories
    i = 0 # Entry index in the DataFrame
    prepep_sweep_list = []
    postpep_sweep_list = []
    for d in subdirs:
        subdir_path = os.path.join(bin_dir, d)
        file_list = os.listdir(subdir_path)
        file_list.sort()

        # Iterate over all files in the subdirectory
        for f in file_list:
            # Sanity check file date and run metadata against spreadsheet
            fdate = datetime.strptime("/".join(f.split(".")[0].rsplit("_", 3)[1:]),
                                      "%m/%d/%Y").date()
            frun = int(f.split("_", 1)[0][1])
            sheet_date = datetime.strptime(sheet_df["date"][i], "%m/%d/%Y").date()
            sheet_run = sheet_df["run_num"][i]
            if fdate != sheet_date or frun != sheet_run:
                print(f"""Entry in spreadsheet does not match. Sheet has \
{sheet_date} {sheet_run} but next file has {fdate} {frun}""")
                sys.exit(1)

            # Check experiment type, skip if not matching
            key = _derive_experiment_key(sheet_df, i)
            if key != type_key:
                i += 1
                continue

            sweep_idx = _find_first_transient_cell(sheet_df, i)

            # Skip confounded experiment
            if sheet_df["data_OK"][i] == False:
                if verbose:
                    run_id = f"{sheet_df["date"][i]} {sheet_df["run_num"][i]}"
                    print("Skipping confounded experiment " + run_id)
            else:
                fpath = os.path.join(subdir_path, f)
                data = binarywave.load(fpath)['wave']['wData']

                # Iterate over all sweeps in the file
                for d in range(data.shape[1]):
                    df_sweep_name = "transient" + ("+" if sweep_idx >= 0 else "") + str(sweep_idx)
                    sweep_idx += 1
                    # Skip if the sweep is marked as nan in the sheet
                    if sheet_df[df_sweep_name][i] == np.nan:
                        continue

                    if sweep_idx < 0:
                        prepep_sweep_list.append(data[:, d])
                    else:
                        postpep_sweep_list.append(data[:, d])

            i += 1 # Increment DataFrame index

    prepep_array = np.array(prepep_sweep_list)
    postpep_array = np.array(postpep_sweep_list)
    return prepep_array, postpep_array


def _get_manual_transient_count_single(df : pd.DataFrame, idx : int):
    """Get pre-PEP and post-PEP transient counts for each sweep in a given
        experiment, as tallied in the master sheet.
    
    Args:
        sheet_df (DataFrame): The DataFrame containing metadata from the master sheet.
        idx (int): Index of the entry in sheet_df.
    
    Returns:
        tuple:
            expt_type (str): string key of the experiment's type.
            prepep_transients (list): list of transient counts for pre-PEP sweeps.
            postpep_transients (list): list of transient counts for post-PEP sweeps.
    """
    expt_type = _derive_experiment_key(df, idx)
    sweep_idx = _find_first_transient_cell(df, idx)
    
    sweep_name = "transient"+("" if sweep_idx < 0 else "+")+str(sweep_idx)
    prepep_transients = []
    postpep_transients = []
    while (sweep_idx < MAX_POST_PEP_SWEEP and df[sweep_name][idx] != ""):
        # Make sure nan values don't get added
        if not np.isnan(df[sweep_name][idx]):
            if sweep_idx < 0:
                prepep_transients.append(df[sweep_name][idx])
            else:
                postpep_transients.append(df[sweep_name][idx])
        sweep_idx += 1
        sweep_name = "transient"+("" if sweep_idx < 0 else "+")+str(sweep_idx)

    return expt_type, prepep_transients, postpep_transients

def _init_transient_dict(type_file_path : str):
    """Initialize a dictionary with a list for each experiment type in the
        experiment type list file.

    Args:
        type_file_path (str): File path to load type list from.

    Returns:
        tuple:
            expt_type_list (list[str]): List of experiment types.
            transients_by_expt_type (dict[list]): Dictionary of lists of
                transient counts for pre- and post-PEP segments of each
                experiment type.    
    """
    expt_types = eh.get_expt_types(type_file_path)
    pre_keys, post_keys = eh.get_expt_keys(type_file_path)

    # Pre- vs. post-PEP
    transients_by_expt_type = {k: [] for k in pre_keys} |\
        {k: [] for k in post_keys}

    return expt_types, transients_by_expt_type

def get_manual_transient_count(df : pd.DataFrame, type_list_path : str,
                               aggregate : bool = True):
    """Get pre-PEP and post-PEP transient counts for all experiments, as
        tallied in the master sheet.
    
    Args:
        sheet_df (DataFrame): The DataFrame containing metadata from the master sheet.
        type_list_path (str): File path of experiment type list.
        aggregate (bool): Compute the sum across all sweeps in each run if True,
            otherwise keep per-sweep counts separate.
        
    Returns:
        tuple:
            expt_types (list[str]): list of key strings of experiment types.
            transients (dict[list]): dictionary of lists of transient counts
                for pre- and post-PEP segments of each experiment type.
    """
    row_count = len(df)
    expt_types, transients = _init_transient_dict(type_list_path)
    for i in range(row_count):
        if df["data_OK"][i] == False:
            continue
        expt, prepep, postpep = _get_manual_transient_count_single(df, i)
        if expt not in expt_types:
            continue

        if aggregate:
            prepep = [np.sum(prepep)]
            postpep = [np.sum(postpep)]

        transients[expt + "_pre"].extend(prepep)
        transients[expt + "_post"].extend(postpep)

    return expt_types, transients

def analyze_binaries(sheet_df : pd.DataFrame, bin_dir : str, type_list_path : str, 
                     detect_freq = DRIFT_FILT_FREQ, detect_thresh = STD_THRESH,
                     verbose=False, aggregate : bool = True):
    """Analyze a directory of binary files using a master spreadsheet as metadata.
        Binaries must be organized into subdirectories by date of experiment, and
        subdirectories must be named in the format '[m]m dd yyyy'. Binaries must be
        named starting with 'R[run number]_'.

        Args:
            sheet_df (DataFrame): The processed DataFrame derived from the
                master spreadsheet.
            bin_dir (str): The path to the topmost binary directory.
            type_list_path (str): Path to experiment type list file.
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

    expt_types, transients_by_expt_type = _init_transient_dict(type_list_path)
    # Iterate over all subdirectories
    i = 0 # Entry index in the DataFrame
    for d in subdirs:
        subdir_path = os.path.join(bin_dir, d)
        file_list = os.listdir(subdir_path)
        file_list.sort()

        # Iterate over all files in the subdirectory
        for f in file_list:
            # Sanity check file date and run metadata against spreadsheet
            fdate = datetime.strptime("/".join(f.split(".")[0].rsplit("_", 3)[1:]),
                                      "%m/%d/%Y").date()
            frun = int(f.split("_", 1)[0][1])
            sheet_date = sheet_df["date"][i]
            sheet_run = sheet_df["run_num"][i]
            if fdate != sheet_date or frun != sheet_run:
                print(f"""Entry in spreadsheet does not match. Sheet has \
{sheet_date} {sheet_run} but next file has {fdate} {frun}""")
                sys.exit(1)

            start_sweep_count = _find_first_transient_cell(sheet_df, i)

            if sheet_df["data_OK"][i] == False:
                if verbose:
                    run_id = f"{sheet_df["date"][i]} {sheet_df["run_num"][i]}"
                    print("Skipping confounded experiment " + run_id)
            else:
                fpath = os.path.join(subdir_path, f)
                pre_counted, post_counted = get_binary_transients(fpath, df=sheet_df, entry_idx=i,
                                                                    df_start_pos=start_sweep_count,
                                                                    detect_freq=detect_freq,
                                                                    detect_thresh=detect_thresh,
                                                                    verbose=verbose, show=False)
                
                if aggregate:
                    pre_counted = [np.sum(pre_counted)]
                    post_counted = [np.sum(post_counted)]

                # Add transient counts to appropriate lists
                # Skip if key not in bins dict
                key = _derive_experiment_key(sheet_df, i)
                if key in expt_types:
                    transients_by_expt_type[key + "_pre"].extend(pre_counted)
                    transients_by_expt_type[key + "_post"].extend(post_counted)
                #print(len(pre_counted), len(post_counted))

                if verbose:
                    print("Expected transients:", sheet_df["total_transients"][i])
                    print(f"Pre-PEP transients: {int(pre_counted[0])}")
                    print(f"Post-PEP transients: {int(post_counted[0])}")
                
            i += 1 # Increment DataFrame index
            if verbose:
                print()

    return expt_types, transients_by_expt_type

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--fname', '-f', type=str, default=None,
                        help='Path to ibw binary file.')
    parser.add_argument('--in-dir', '-i', type=str, default=DEFAULT_IN_DIR,
                        help="""Path to directory of ibw binary files. 
Subdirectories must have naming scheme [m]m dd yyyy.""")
    parser.add_argument('--sweep', default=0, type=int,
                        help='The 1-indexed sweep number to examine. Leave out for all sweeps')
    parser.add_argument('--sheet-path', default=DEFAULT_SHEET_PATH, type=str,
                        help="""Path to master spreadsheet file to supply ancillary data for binary analysis.""")
    parser.add_argument('--show', default=False, action='store_true',
                        help='Show plots after completion.')
    parser.add_argument('--out-dir', '-o', default=DEFAULT_OUT_DIR, type=str,
                        help='Directory to save figures to.')
    parser.add_argument('-v', '--verbose', action='count',
                        default=0, help='Verbosity of operations. Type vv for level 2 verbosity.')
    parser.add_argument('--mode', '-m', type=str, default=None,
                        required=True,
                        help="""Plotting mode.
    \nplot_sweep: Plot LFPs from a single file specified by -f or --fname.
    Pass --sweep to specify one sweep, otherwise all sweeps will be plotted. 
    \nscatter_col: Plot scatter column plot of pre- and post-PEP spike counts by
    experiment type.
    \nmean_cv: Plot two-row bar graph of mean and CV of spike counts by experiment
    type, plotting both pre- and post-PEP.
    \nauto_v_man: Plot a scatter plot of the automatic spike tally against the
    manual one.
    \npower: Plot mean pre- and post-PEP power spectral density for a given run.
    \npower_by_type: Plot mean pre- and post-PEP power spectral density of all runs of
    a given type.
    \nmean_spikes: Print mean pre- and post-PEP spike counts for each experiment type.
    \nwilcoxon: Print statistic value and p-value from Wilcoxon signed-rank test for
    each experiment type.
    \nchange_score: Print pre-to-post change scores for each experiment of each type.
    \n""")
    parser.add_argument('--type', '-t', type=str, default=None,
                        help="Type of experiment for analysis. Used by power_by_type.")
    parser.add_argument('--man', action='store_true', default=False,
                        help="Use manual counts from master sheet.")
    parser.add_argument('--expt-list', '-e', type=str, default=DEFAULT_EXPERIMENT_FILE)

    args = parser.parse_args()
    fname = args.fname
    in_dir = args.in_dir
    sheet_path = args.sheet_path
    out_dir = args.out_dir
    sweep = args.sweep - 1
    show = args.show
    verbose = args.verbose
    mode = args.mode
    etype = args.type
    manual = args.man
    type_list_path = args.expt_list
    sheet_df = sheet_parser.parse_sheet(sheet_path)

    if mode == None:
        print("Please specify a mode with --mode.")
        sys.exit(1)

    if mode != 'detect_param_sweep' and mode != 'specgram_pdf':
        if mode == 'auto_v_man':
            expt_types, transients_auto = analyze_binaries(sheet_df, in_dir, type_list_path,
                                                            verbose=(verbose > 1),
                                                            aggregate=False)
            expt_types, transients_man = get_manual_transient_count(sheet_df, type_list_path,
                                                                    aggregate=False)
        else:
            if manual:
                expt_types, transients = get_manual_transient_count(sheet_df, type_list_path)
                title_suffix = "(Manual Counting)"
                out_suffix = "_man"
            else:
                expt_types, transients = analyze_binaries(sheet_df, in_dir, type_list_path,
                                                        verbose=(verbose > 1))
                title_suffix = "(Automated Counting)"
                out_suffix = "_auto"
            expt_types, split_point = eh.separate_prepost_columns(expt_types, transients)
    
    match mode:
        case 'count':
            if fname == None:
                print("Please specify a file to plot with --fname")
                sys.exit(1)
            else:
                start_sweep, idx = _find_first_transient_cell_by_file(sheet_df, fname)
                prepep, postpep = get_binary_transients(fname, df=sheet_df, entry_idx=idx,
                                                        verbose=(verbose > 0),
                                                        df_start_pos=start_sweep, scan_for_nan=False)
                print("Pre-PEP transients:", prepep)
                print("Post-PEP transients:", postpep)
        case 'plot_sweep':
            if fname == None:
                print("Please specify a file to plot with --fname")
                sys.exit(1)
            else:
                get_binary_transients(fname, disp_sweep=sweep, show=show,
                                       verbose=(verbose > 0), out_dir=out_dir)
        case 'scatter_col':
            ephys_plots.plot_scatter_columns(expt_types, split_point,
                                                transients, show=show, out_dir=out_dir,
                                                out_suffix=out_suffix,
                                                title="Spikes by Experiment " + title_suffix)
        case 'mean_cv':
            ephys_plots.plot_mean_cv_bar(expt_types, split_point,
                                            transients, out_dir=out_dir,
                                            show=show, title="Spike Stats " + title_suffix,
                                            out_suffix=out_suffix)
        case 'auto_v_man':
            ephys_plots.plot_auto_v_man(transients_man, transients_auto,
                                        type_list_path=type_list_path,
                                        out_dir=out_dir, show=show)
        case 'power':
            if fname == None:
                print("Please specify a file to plot with --fname")
                sys.exit(1)
            else:
                analysis_plots.plot_ephys_power_spec(out_dir, [binarywave.load(fname)['wave']['wData'][:, sweep]])
        case 'power_by_type':
            if etype == None:
                print("Please specify a type of experiment to analyze with --type or -t")
                sys.exit(1)
            else:
                prepep_array, postpep_array = get_sweeps_by_expt_type(sheet_df, in_dir, etype, verbose)
                if prepep_array.size == 0 or postpep_array.size == 0:
                    print("Specified type had 0 matches. Make sure you spelled it correctly.")
                    sys.exit(1)
                analysis_plots.plot_ephys_mean_power_spec(out_dir, [prepep_array, postpep_array], fmax=40,
                                                          fname="ephys_mean_power_"+etype,
                                                          labels=["Pre-PEP", "Post-PEP"])
        case 'specgram_pdf':
            if fname == None: 
                print("Please specify a file to plot with --fname")
                sys.exit(1)
            else:
                start_sweep_count, idx = _find_first_transient_cell_by_file(sheet_df, fname)
                lfp_list = get_lfp_list(fname, df=sheet_df, entry_idx=idx,
                                        df_start_pos=start_sweep_count)
                out_name = os.path.join(out_dir, os.path.splitext(os.path.basename(fname))[0]+\
                                        "_sweep_spectrogram.pdf")
                ephys_plots.ephys_spectrogram_pdf(out_name, lfp_list, EPHYS_FS, fmax=100)
        case 'mean_spikes':
            # Print mean transient counts
            print("Mean spike counts " + title_suffix)
            print(f"{"Type":<20}{"Pre":>20}{"Post":>20}")
            for expt in expt_types:
                pre_k = expt+"_pre"
                post_k = expt+"_post"
                pre_mean = np.mean(transients[pre_k])
                post_mean = np.mean(transients[post_k])
                print(f"{expt:<20}{f"{pre_mean:.4f}":>20}{f"{post_mean:.4f}":>20}")
        case 'wilcoxon':
            print("Wilcoxon tests " + title_suffix)
            print(f"{"Type":<20}{"Statistic":>20}{"p-value":>20}")
            stat_list, pvalue_list = eh.classwise_wilcoxon(transients, expt_types)
            for expt, statistic, pvalue in zip(expt_types, stat_list, pvalue_list):
                print(f"{expt:<20}{f"{statistic:.4f}":>20}{f"{pvalue:.4f}":>20}")
        case 'change_score':
            # Print change scores for each experiment in each type
            print("Change scores " + title_suffix)
            for expt in expt_types:
                pre_k = expt+"_pre"
                post_k = expt+"_post"
                print(expt)
                for pre_ent, post_ent in zip(transients[pre_k], transients[post_k], strict=True):
                    print(f"{(post_ent+1) / (pre_ent+1) - 1:.4f} ", end='')
                print()
        case 'detect_param_sweep':
            param_step_count = 11
            detect_freq_start = 75
            detect_freq_end = 225
            detect_freq_step = (detect_freq_end - detect_freq_start) / (param_step_count - 1)
            detect_thresh_start = 6
            detect_thresh_end = 10
            detect_thresh_step = (detect_thresh_end - detect_thresh_start) / (param_step_count - 1)
            rsquare_results_pre = np.zeros((param_step_count, param_step_count))
            rsquare_results_post = np.zeros((param_step_count, param_step_count))
            freq_range = np.arange(detect_freq_start, detect_freq_end, detect_freq_step)
            thresh_range = np.arange(detect_thresh_start, detect_thresh_end, detect_thresh_step)
            _, man_transients = get_manual_transient_count(sheet_df, type_list_path, aggregate=False)
            for i, freq in enumerate(freq_range):
                for j, thresh in enumerate(thresh_range):
                    _, auto_transients = analyze_binaries(sheet_df, in_dir, type_list_path,
                                                     detect_freq=freq, detect_thresh=thresh,
                                                     verbose=(verbose > 1),
                                                     aggregate=False)
                    auto_pre, auto_post = eh.split_prepost_transients(auto_transients)
                    man_pre, man_post = eh.split_prepost_transients(man_transients)
                    _, _, r_pre, _, _ = linregress(man_pre, auto_pre)
                    _, _, r_post, _, _ = linregress(man_post, auto_post)
                    rsquare_results_pre[i, j] = r_pre**2
                    rsquare_results_post[i, j] = r_post**2
            print("Pre-PEP R squared results")
            print(rsquare_results_pre)
            max_ind = np.argmax(rsquare_results_pre)
            print("Max R squared:", np.max(rsquare_results_pre),
                  "from parameters freq", freq_range[max_ind // param_step_count],
                  "thresh", thresh_range[max_ind % param_step_count])
            print("Post-PEP R squared results")
            print(rsquare_results_post)
            max_ind = np.argmax(rsquare_results_post)
            print("Max R squared:", np.max(rsquare_results_post),
                  "from parameters freq", freq_range[max_ind // param_step_count],
                  "thresh", thresh_range[max_ind % param_step_count])
        case _:
            print("Mode not recognized. Type --help for mode list.")
               
if __name__ == "__main__":
    main()