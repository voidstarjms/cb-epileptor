import numpy as np
from matplotlib import pyplot as plt
from matplotlib import ticker
import os
from scipy.stats import linregress, wilcoxon

def plot_lfp(x : np.array, y_raw : np.array, y_filtered : np.array, y_processed : np.array,
             std_thresh : float, infile : str, disp_sweep : int, out_dir : str = None):
    """Plot the flattened unprocessed LFP and fully processed LFP on the same x-axis.

    Args:
        x (np.array): Time axis.
        y_raw (np.array): Flattened LFP reading (baseline drift removed).
        y_processed (np.array): Fully processed LFP reading (transients only).
        infile (str): Input file name. Supplies date in the figure title.
        disp_sweep (int): Which sweep from the experiment to plot. -1 is all.
        out_dir (str): Path to save the figures to. Pass None to forgo saving.
    
    Returns: void
    """
    # Plot all LFPs overlaid
    if disp_sweep == -1:
        fig = plt.figure()
        ax = fig.add_axes((0.1, 0.1, 0.8, 0.8))
        fig.set_size_inches(10, 6)
        infile_parts = infile.split(os.sep)[-1].split(".")[0].split("_")
        sweep_num = f"Sweep {disp_sweep+1}" if disp_sweep >= 0 else "All Sweeps"
        fdate = f"{infile_parts[-3]}-{infile_parts[-2]}-{infile_parts[-1]} {infile_parts[0]} {sweep_num}"
        fig.suptitle(f"{fdate}")
        for i in range(y_raw.shape[1]):
            ax.plot(x, y_raw[:, i], alpha=0.4)

        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Raw LFP (V)")
    # PLot a single LFP and the transients detected from it
    else:
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
        fig.set_size_inches(18, 10)
        infile_parts = infile.split(os.sep)[-1].split(".")[0].split("_")
        sweep_num = f"Sweep {disp_sweep+1}" if disp_sweep >= 0 else "All Sweeps"
        fdate = f"{infile_parts[-3]}-{infile_parts[-2]}-{infile_parts[-1]} {infile_parts[0]} {sweep_num}"
        fig.suptitle(f"{fdate}", fontsize=22)
        ax1.plot(x * 1000, y_raw[:, disp_sweep])
        ax2.plot(x * 1000, y_filtered[:, disp_sweep])
        ax3.plot(x * 1000, y_processed[:, disp_sweep])

        ax3.xaxis.set_major_locator(ticker.MaxNLocator(nbins=21))

        ax1.set_ylabel("Raw LFP (V)", fontsize=14)
        ax2.set_ylabel("Filtered\nLFP (V)", fontsize=14)
        ax3.set_ylabel("Detected\nTransients (V)", fontsize=14)
        ax3.set_xlabel("Time (ms)", fontsize=14)
        
        ax1.tick_params(axis='both', which='major', labelsize=12)
        ax2.tick_params(axis='both', which='major', labelsize=12)
        ax3.tick_params(axis='both', which='major', labelsize=12)

        ax1.yaxis.get_offset_text().set_fontsize(12)
        ax2.yaxis.get_offset_text().set_fontsize(12)
        ax3.yaxis.get_offset_text().set_fontsize(12)

        ax2.axhline(std_thresh[disp_sweep], c='r')
        ax3.axhline(std_thresh[disp_sweep], c='r', label="Threshold")

        plt.legend(prop={'size': 12})

    #plt.xlabel(w['dimension_units'].decode().strip('\x00'))
    #plt.ylabel(w['data_units'].decode().strip('\x00'))
    if out_dir != None:
        outfile = "ephys_lfp_" + ("all_sweeps" if disp_sweep == -1 else str(disp_sweep+1))
        plt.savefig(os.path.join(out_dir, outfile))

    plt.show()

def _filter_entryless_columns(expt_types : list, transients : dict):
    """Delete experiment types with empty transient count lists.

    Args:
        expt_types (list[str]): List of experiment type names, keys to the transients dict.
        transients (dict[list]): Dictionary of lists containing transient counts for each
            experiment of each experiment type.
        
    Returns:
        tuple (prepep_keys, postpep_keys, col_ids, split_point)
            prepep_keys: Keys corresponding to pre-PEP transient counts.
            postpep_keys: Keys corresponding to post-PEP transient counts.
            col_ids: List of column positions.
            split_point: The column at which the second arm of the experiment type set
                starts in the column list (This will usually be WT vs. BTBR or vehicle vs.
                CBD).
    """
    # Check for columns with no entries
    no_entry_keys = []
    max_col_num = len(expt_types)
    col_num = max_col_num
    split_point = col_num / 2
    for i, k in enumerate(transients.keys()):
        if len(transients[k]) == 0:
            # Add key to list of keys with no entries
            no_entry_keys.append(k)
            # Decrease the number of columns
            if i <= col_num:
                col_num -= 1
            # Keep track of the first BTBR column
            # (WTF why is the +2 needed?)
            if i < split_point + 2:
                split_point -= 1
    
    # Delete all entryless columns
    for k in no_entry_keys:
        # Delete the key from the transients dict
        del transients[k]

    # Delete the corresponding experiment name from the name list
    # (stop at col_num to ignore distinct keys for pre and post)
    new_expt_types = []
    for e in expt_types:
        if e+"_pre" not in no_entry_keys:
            new_expt_types.append(e)
    expt_types = new_expt_types

    # Generate column positions
    col_ids = np.arange(col_num)

    # Construct lists of keys, scatter x positions, and scatter y positions for prepep and postpep
    prepep_keys = list(transients.keys())[:col_num]
    postpep_keys = list(transients.keys())[col_num:]

    return (prepep_keys, postpep_keys, expt_types, col_ids, split_point)

def _make_column_labels(expt_types : list):
    # Generate column labels
    col_labels = []
    for t in expt_types:
        label = ""
        label += "BTBR" if "btbr" in t else "WT"
        if "control" in t:
            label += "\nControl"
        if "cbd" in t:
            label += "\nCBD"
            if "10uM" in t:
                label += " (10uM)"
            else:
                label += " (100uM)"
        if "picro" in t:
            label += "\nPicrotoxin"
        if "bc" in t:
            label += "\nBC"
        col_labels.append(label)

    return col_labels

def plot_mean_cv_bar(expt_types : list, transients : dict, out_dir : str = None, show : bool = False,
                     title : str = "", out_suffix : str = ""):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 12), sharex=True)

    prepep_keys, postpep_keys, col_ids, split_point = _filter_entryless_columns(
        expt_types, transients)
    col_num = len(col_ids)

    offset = 0.2
    prepep_arrays = [np.array(transients[k]) for k in prepep_keys]
    postpep_arrays = [np.array(transients[k]) for k in postpep_keys]
    prepep_means = [np.mean(a) for a in prepep_arrays]
    postpep_means = [np.mean(a) for a in postpep_arrays]
    prepep_cvs = [np.std(a) / prepep_means[i] if prepep_means[i] != 0 else 0
                  for i, a in enumerate(prepep_arrays)]
    postpep_cvs = [np.std(a) / postpep_means[i] if prepep_means[i] != 0 else 0
                   for i, a in enumerate(postpep_arrays)]
    ax1.bar(col_ids - offset, prepep_means, width=0.3)
    ax1.bar(col_ids + offset, postpep_means, width=0.3)
    ax1.set_ylabel("Mean Transient Count", fontsize=20)
    ax1.tick_params(axis='both', which='major', labelsize=15)
    for i in range(0, col_num - 1):
        ax1.axvline(i + 0.5, c='k', lw=1)
    if split_point > 0:
        ax1.axvline(split_point - 0.5, c='k', lw=3)
    ax1.grid(visible=True, axis='y', which='major')

    ax2.bar(col_ids - offset, prepep_cvs, label="Pre-PEP", width=0.3)
    ax2.bar(col_ids + offset, postpep_cvs, label="Post-PEP", width=0.3)
    ax2.set_ylabel("CV of Transient Count", fontsize=20)
    ax2.set_xlabel("Experiment Type", fontsize=20)
    ax2.tick_params(axis='both', which='major', labelsize=15)
    for i in range(0, col_num - 1):
        ax2.axvline(i + 0.5, c='k', lw=1)
    if split_point > 0:
        ax2.axvline(split_point - 0.5, c='k', lw=3)
    ax2.grid(visible=True, axis='y', which='major')

    plt.xticks(col_ids, _make_column_labels(expt_types))
    plt.legend()
    fig.suptitle(title, fontsize=30)

    plt.savefig(os.path.join(out_dir, "ephys_mean_cv"+out_suffix))
    if show:
        plt.show()

def plot_scatter_columns(expt_types : list, transients : dict, out_dir : str = None, show : bool = False, out_suffix : str = "",
                         title : str = None):
    """Scatter plot organized by column, with each column representing a class of experiment.
        Pre-PEP transient counts are to the left of the column, post-PEP to the right, with
        lines connecting the pre and post values for a given experiment.

    Args:
        expt_types (list[str]): List of experiment type names, keys to the transients dict.
        transients (dict[list]): Dictionary of lists containing transient counts for each
            experiment of each experiment type.
        out_dir (str): Directory to save figures to. Pass none to forgo saving.
        show (bool): Show the plot.
        out_suffix (str): Suffix for the file name, appended before the extension.

    Returns: void
    """

    fig = plt.figure(figsize=(20, 10))
    fig.add_axes((0.075, 0.125, 0.85, 0.75))
    ax1 = fig.axes[0]

    prepep_keys, postpep_keys, expt_types, col_ids, split_point = _filter_entryless_columns(
        expt_types, transients)
    col_num = len(col_ids)

    jitter = 0.05
    offset = 0.25

    # Generate dots for pre- and post-PEP data. Jitter x positions for ease of viewing
    prepep_x_list = [col + np.random.uniform(-jitter, jitter) - offset
                     for col, k in enumerate(prepep_keys)
                     for _ in transients[k]]
    prepep_y_list = [v for k in prepep_keys for v in transients[k]]
    postpep_x_list = [col + np.random.uniform(-jitter, jitter) + offset
                      for col, k in enumerate(postpep_keys)
                      for _ in transients[k]]
    postpep_y_list = [v for k in postpep_keys for v in transients[k]]

    # Generate column labels
    col_labels = _make_column_labels(expt_types)

    # Plot column scatter plot of transient counts
    ax1.scatter(prepep_x_list, prepep_y_list, color=(0.0, 0.0, 1.0, 0.5),
                marker='o', s=80, label="Pre-PEP")
    ax1.set_ylabel("# of Transients", fontsize=22)
    ax1.tick_params(axis='both', which='major', labelsize=15)
    for i in range(0, col_num - 1):
        ax1.axvline(i + 0.5, c='k', lw=1)
    if split_point > 0:
        ax1.axvline(split_point - 0.5, c='k', lw=3)
    ax1.grid(visible=True, axis='y', which='major')

    ax1.scatter(postpep_x_list, postpep_y_list, color=(1.0, 0.0, 0.0, 0.5),
                marker='o', s=80, label="Post-PEP")
    for pre_x, pre_y, post_x, post_y in zip(prepep_x_list, prepep_y_list, postpep_x_list, postpep_y_list):
        ax1.plot([pre_x, post_x], [pre_y, post_y], c='k')
    
    # Column labels
    ax1.set_xlabel("Experiment Type", fontsize=22)
    plt.xticks(col_ids, col_labels, fontsize=15)
    ax1.legend(prop={'size': 12})
    fig.suptitle(title, fontsize=30)

    # Print mean transient counts
    print("Mean transient counts")
    print(f"{"Type":<20}{"Pre":>20}{"Post":>20}")
    for pre_k, post_k in zip(prepep_keys, postpep_keys):
        pre_mean = np.mean(transients[pre_k])
        post_mean = np.mean(transients[post_k])
        print(f"{pre_k.rsplit("_", 1)[0]:<20}{f"{pre_mean:.4f}":>20}{f"{post_mean:.4f}":>20}")
    print()

    # Compute Wilcoxon signed-rank test for each column
    print("Wilcoxon tests")
    print(f"{"Type":<20}{"Statistic":>20}{"p-value":>20}")
    pvalue_text_y = 1.02
    pvalue_num_y = ax1.get_ylim()[1] * pvalue_text_y
    pvalue_text_col = 'b'
    plt.text(-0.075, pvalue_text_y, "p-value", transform=ax1.transAxes,
             fontsize=20, c=pvalue_text_col)
    for i, pre_k, post_k in zip(ax1.get_xticks(), prepep_keys, postpep_keys):
        x = np.array(transients[pre_k])
        y = np.array(transients[post_k])
        wilcoxon_results = wilcoxon(y - x)
        statistic = wilcoxon_results.statistic
        pvalue = wilcoxon_results.pvalue
        plt.text(i, pvalue_num_y, f"{pvalue:.4f}",
                 transform=ax1.transData, horizontalalignment='center', fontsize=15,
                 c=pvalue_text_col)
        print(f"{pre_k.rsplit("_", 1)[0]:<20}{f"{statistic:.4f}":>20}{f"{pvalue:.4f}":>20}")

    if out_dir != None:
        plt.savefig(os.path.join(out_dir, "ephys_scatter_columns"+out_suffix))

    if show:
        plt.show()

def plot_auto_v_man(expt_types : list, man_transients : dict, auto_transients : dict, out_dir : str = None,
                    show : bool = False):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 9))

    man_prepep_keys, man_postpep_keys, _, _, _ = _filter_entryless_columns(expt_types, man_transients)
    auto_prepep_keys, auto_postpep_keys, _, _, _ = _filter_entryless_columns(expt_types, auto_transients)

    man_pre_transients = [e for k in man_prepep_keys for e in man_transients[k]]
    auto_pre_transients = [e for k in auto_prepep_keys for e in auto_transients[k]]
    man_post_transients = [e for k in man_postpep_keys for e in man_transients[k]]
    auto_post_transients = [e for k in auto_postpep_keys for e in auto_transients[k]]

    # Plot pre-PEP comparison
    print(len(man_pre_transients), len(auto_pre_transients))
    print(len(man_post_transients), len(auto_post_transients))
    ax1.scatter(man_pre_transients, auto_pre_transients)
    ax1.tick_params(axis='both', which='major', labelsize=15)
    ax1.set_xlabel("# Transients (Manual)", fontsize=20)
    ax1.set_ylabel("# Transients (Automated)", fontsize=20)
    ax1.set_title("Pre-PEP", fontsize=20)
    upper_lim = max(ax1.get_xlim()[1], ax1.get_ylim()[1])
    ax1.set_xlim(right=upper_lim)
    ax1.set_ylim(top=upper_lim)
    ax1.set_aspect('equal', adjustable='box')
    ax1.plot([0, 1], [0, 1], transform=ax1.transAxes, c='k')

    # Plot post-PEP comparison
    ax2.scatter(man_post_transients, auto_post_transients)
    ax2.tick_params(axis='both', which='major', labelsize=15)
    ax2.set_xlabel("# Transients (Manual)", fontsize=20)
    ax2.set_ylabel("# Transients (Automated)", fontsize=20)
    ax2.set_title("Post-PEP", fontsize=20)
    upper_lim = max(ax2.get_xlim()[1], ax2.get_ylim()[1])
    ax2.set_xlim(right=upper_lim)
    ax2.set_ylim(top=upper_lim)
    ax2.set_aspect('equal', adjustable='box')
    ax2.plot([0, 1], [0, 1], transform=ax2.transAxes, c='k')

    fig.suptitle("Automated v. Manual Transient Counting", fontsize=30)
    fig.tight_layout(pad=2)

    _, _, r_value, _, _ = linregress(man_pre_transients, auto_pre_transients)
    print(f"Pre-PEP r squared: {r_value**2}")
    _, _, r_value, _, _ = linregress(man_post_transients, auto_post_transients)
    print(f"Post-PEP r squared: {r_value**2}")

    if out_dir:
        plt.savefig(os.path.join(out_dir, "ephys_auto_man_comparison"))
    if show:
        plt.show()