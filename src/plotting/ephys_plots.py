import numpy as np
from matplotlib import pyplot as plt
import os
from scipy.stats import linregress, wilcoxon

def plot_lfp(x : np.array, y_raw : np.array, y_processed : np.array,
             infile : str, disp_sweep : int, out_dir : str = None):
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
        fig.set_size_inches(7, 6)
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
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        fig.set_size_inches(14, 6)
        infile_parts = infile.split(os.sep)[-1].split(".")[0].split("_")
        sweep_num = f"Sweep {disp_sweep+1}" if disp_sweep >= 0 else "All Sweeps"
        fdate = f"{infile_parts[-3]}-{infile_parts[-2]}-{infile_parts[-1]} {infile_parts[0]} {sweep_num}"
        fig.suptitle(f"{fdate}")
        ax1.plot(x, y_raw[:, disp_sweep])
        ax2.plot(x, y_processed[:, disp_sweep])

        ax2.set_xlabel("Time (s)")
        ax1.set_ylabel("Raw LFP (V)")
        ax2.set_ylabel("Processed LFP\n(V relative to mean)")

    #plt.xlabel(w['dimension_units'].decode().strip('\x00'))
    #plt.ylabel(w['data_units'].decode().strip('\x00'))
    if out_dir != None:
        outfile = "ephys_lfp_" + "all_sweeps" if disp_sweep == -1 else str(disp_sweep)
        plt.savefig(os.path.join(out_dir, outfile))

    plt.show()


def plot_scatter_columns(expt_types : list, transients : dict, out_dir : str = None, show : bool = False):
    """Scatter plot organized by column, with each column representing a class of experiment.
        Pre-PEP transient counts are to the left of the column, post-PEP to the right, with
        lines connecting the pre and post values for a given experiment.

    Args:
        expt_types (list[str]): List of experiment type names, keys to the transients dict.
        transients (dict[list]): Dictionary of lists containing transient counts for each
            experiment of each experiment type.
        out_dir (str): Directory to save figures to. Pass none to forgo saving.
        show (bool): Show the plot.

    Returns: void
    """

    fig = plt.figure(figsize=(20, 12))
    fig.add_axes((0.05, 0.1, 0.9, 0.8))
    ax1 = fig.axes[0]

    # Check for columns with no entries
    no_entry_keys = []
    max_col_num = len(expt_types)
    col_num = max_col_num
    btbr_start_col = col_num / 2
    for i, k in enumerate(transients.keys()):
        if len(transients[k]) == 0:
            # Add key to list of keys with no entries
            no_entry_keys.append((i, k))
            # Shrink the starting position of the prepep region as it decreases in size
            if i < col_num:
                col_num -= 1
    
    # Delete all entryless columns
    for i, k in no_entry_keys:
        # Delete the key from the transients dict
        del transients[k]
        # Delete the corresponding experiment name from the name list
        # (stop at max_col_num to ignore distinct keys for pre and post)
        if i < max_col_num:
            del expt_types[i]
        # Keep track of the first BTBR column
        if i < btbr_start_col:
            btbr_start_col -= 1

    # Generate column positions
    col_ids = np.arange(col_num)

    # Construct lists of keys, scatter x positions, and scatter y positions for prepep and postpep
    prepep_keys = list(transients.keys())[:col_num]
    postpep_keys = list(transients.keys())[col_num:]
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

    # Compute Wilcoxon signed-rank test for each column
    print("Wilcoxon tests")
    print(f"{"Type":<15}{"Statistic":>15}{"p-value":>15}")
    for pre_k, post_k in zip(prepep_keys, postpep_keys):
        x = np.array(transients[pre_k])
        y = np.array(transients[post_k])
        wilcoxon_results = wilcoxon(y - x)
        statistic = wilcoxon_results.statistic
        pvalue = wilcoxon_results.pvalue
        print(f"{pre_k.rsplit("_", 1)[0]:<15}{f"{statistic:.4f}":>15}{f"{pvalue:.4f}":>15}")

    # Generate column labels
    col_labels = []
    for t in expt_types:
        label = ""
        label += "BTBR" if "btbr" in t else "WT"
        if "control" in t:
            label += "\nControl"
        if "cbd" in t:
            label += "\nCBD"
        if "picro" in t:
            label += "\nPicrotoxin"
        if "bc" in t:
            label += "\nBC"
        col_labels.append(label)

    # Plot column scatter plot of transient counts
    ax1.scatter(prepep_x_list, prepep_y_list, color=(0.0, 0.0, 1.0, 0.5),
                marker='o', s=80, label="Pre-PEP")
    ax1.set_ylabel("# of Transients", fontsize=20)
    ax1.tick_params(axis='both', which='major', labelsize=15)
    for i in range(0, col_num - 1):
        ax1.axvline(i + 0.5, c='k', lw=1)
    if btbr_start_col > 0:
        ax1.axvline(btbr_start_col - 0.5, c='k', lw=3)
    ax1.grid(visible=True, axis='y', which='major')

    ax1.scatter(postpep_x_list, postpep_y_list, color=(1.0, 0.0, 0.0, 0.5),
                marker='o', s=80, label="Post-PEP")
    for pre_x, pre_y, post_x, post_y in zip(prepep_x_list, prepep_y_list, postpep_x_list, postpep_y_list):
        ax1.plot([pre_x, post_x], [pre_y, post_y], c='k')
    
    # Column labels
    ax1.set_xlabel("Experiment Type", fontsize=20)
    plt.xticks(col_ids, col_labels)
    ax1.legend()

    if out_dir != None:
        plt.savefig(os.path.join(out_dir, "ephys_scatter_columns"))

    if show:
        plt.show()