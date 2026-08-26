import numpy as np
from matplotlib import pyplot as plt
from matplotlib import ticker
import os
from scipy.stats import linregress, wilcoxon, f_oneway
import sys
sys.path.append("..")
from ephys import ephys_util as eh

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
        ax3.set_ylabel("Detected\nSpikes (V)", fontsize=14)
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

def plot_mean_cv_bar(expt_types : list, split_point : int,
                     transients : dict, out_dir : str = None, show : bool = False,
                     title : str = "", out_suffix : str = ""):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 12), sharex=True)

    prepep_keys = [ent+"_pre" for ent in expt_types]
    postpep_keys = [ent+"_post" for ent in expt_types]
    col_num = len(expt_types)
    col_ids = np.arange(col_num)

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
    ax1.set_ylabel("Mean Spike Count", fontsize=20)
    ax1.tick_params(axis='both', which='major', labelsize=15)
    for i in range(0, col_num - 1):
        ax1.axvline(i + 0.5, c='k', lw=1)
    if split_point > 0:
        ax1.axvline(split_point - 0.5, c='k', lw=3)
    ax1.grid(visible=True, axis='y', which='major')

    ax2.bar(col_ids - offset, prepep_cvs, label="Pre-PEP", width=0.3)
    ax2.bar(col_ids + offset, postpep_cvs, label="Post-PEP", width=0.3)
    ax2.set_ylabel("CV of Spike Count", fontsize=20)
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

def plot_scatter_columns(expt_types : list, split_point : int,
                         transients : dict, out_dir : str = None, show : bool = False,
                         out_suffix : str = "", title : str = None):
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

    prepep_keys = [ent+"_pre" for ent in expt_types]
    postpep_keys = [ent+"_post" for ent in expt_types]
    col_num = len(expt_types)
    col_ids = np.arange(col_num)

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
    ax1.set_ylabel("# of Spikes", fontsize=22)
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

    # Compute Wilcoxon signed-rank test for each column
    _, pvalue_list = eh.classwise_wilcoxon(transients, expt_types)
    pvalue_text_y = 1.02
    pvalue_num_y = ax1.get_ylim()[1] * pvalue_text_y
    pvalue_text_col = 'b'
    plt.text(-0.075, pvalue_text_y, "p-value", transform=ax1.transAxes,
             fontsize=20, c=pvalue_text_col)
    for i, pvalue in zip(ax1.get_xticks(), pvalue_list):
        plt.text(i, pvalue_num_y, f"{pvalue:.4f}",
                 transform=ax1.transData, horizontalalignment='center', fontsize=15,
                 c=pvalue_text_col)

    # print()
    # wt_change_scores = np.divide(np.array(transients['wt_cbd_10uM_post'])+1,
    #                              np.array(transients['wt_cbd_10uM_pre'])+1)
    # btbr_change_scores = np.divide(np.array(transients['btbr_cbd_10uM_post'])+1,
    #                                np.array(transients['btbr_cbd_10uM_pre'])+1)
    # fstat, _ = f_oneway(wt_change_scores, btbr_change_scores, nan_policy='omit')
    # print(f"F-stat of WT v. BTBR CBD: {fstat}")

    if out_dir != None:
        plt.savefig(os.path.join(out_dir, "ephys_scatter_columns"+out_suffix))

    if show:
        plt.show()

def plot_auto_v_man(expt_types : list, man_transients : dict, auto_transients : dict, out_dir : str = None,
                    show : bool = False):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 9))

    prepep_keys = [ent+"_pre" for ent in expt_types]
    postpep_keys = [ent+"_post" for ent in expt_types]

    man_pre_transients = [e for k in prepep_keys for e in man_transients[k]]
    auto_pre_transients = [e for k in prepep_keys for e in auto_transients[k]]
    man_post_transients = [e for k in postpep_keys for e in man_transients[k]]
    auto_post_transients = [e for k in postpep_keys for e in auto_transients[k]]

    # Plot pre-PEP comparison
    print(len(man_pre_transients), len(auto_pre_transients))
    print(len(man_post_transients), len(auto_post_transients))
    ax1.scatter(man_pre_transients, auto_pre_transients)
    ax1.tick_params(axis='both', which='major', labelsize=15)
    ax1.set_xlabel("# Spikes (Manual)", fontsize=20)
    ax1.set_ylabel("# Spikes (Automated)", fontsize=20)
    ax1.set_title("Pre-PEP", fontsize=20)
    upper_lim = max(ax1.get_xlim()[1], ax1.get_ylim()[1])
    ax1.set_xlim(right=upper_lim)
    ax1.set_ylim(top=upper_lim)
    ax1.set_aspect('equal', adjustable='box')
    ax1.plot([0, 1], [0, 1], transform=ax1.transAxes, c='k')

    # Plot post-PEP comparison
    ax2.scatter(man_post_transients, auto_post_transients)
    ax2.tick_params(axis='both', which='major', labelsize=15)
    ax2.set_xlabel("# Spikes (Manual)", fontsize=20)
    ax2.set_ylabel("# Spikes (Automated)", fontsize=20)
    ax2.set_title("Post-PEP", fontsize=20)
    upper_lim = max(ax2.get_xlim()[1], ax2.get_ylim()[1])
    ax2.set_xlim(right=upper_lim)
    ax2.set_ylim(top=upper_lim)
    ax2.set_aspect('equal', adjustable='box')
    ax2.plot([0, 1], [0, 1], transform=ax2.transAxes, c='k')

    fig.suptitle("Automated v. Manual Spike Counting", fontsize=30)
    fig.tight_layout(pad=2)

    _, _, r_value, _, _ = linregress(man_pre_transients, auto_pre_transients)
    print(f"Pre-PEP r squared: {r_value**2}")
    _, _, r_value, _, _ = linregress(man_post_transients, auto_post_transients)
    print(f"Post-PEP r squared: {r_value**2}")

    if out_dir:
        plt.savefig(os.path.join(out_dir, "ephys_auto_man_comparison"))
    if show:
        plt.show()