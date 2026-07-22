import numpy as np
from matplotlib import pyplot as plt
import os

def plot_lfp(x : np.array, y_raw : np.array, y_processed : np.array,
             fname : str, disp_sweep : int):
    # Plot all LFPs overlaid
    if disp_sweep == -1:
        fig = plt.figure()
        ax = fig.add_axes((0.1, 0.1, 0.8, 0.8))
        fig.set_size_inches(7, 6)
        fname_parts = fname.split(os.sep)[-1].split(".")[0].split("_")
        sweep_num = f"Sweep {disp_sweep+1}" if disp_sweep >= 0 else "All Sweeps"
        fdate = f"{fname_parts[-3]}-{fname_parts[-2]}-{fname_parts[-1]} {fname_parts[0]} {sweep_num}"
        fig.suptitle(f"{fdate}")
        for i in range(y_raw.shape[1]):
            ax.plot(x, y_raw[:, i], alpha=0.4)

        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Raw LFP (V)")
    # PLot a single LFP and the transients detected from it
    else:
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
        fig.set_size_inches(14, 6)
        fname_parts = fname.split(os.sep)[-1].split(".")[0].split("_")
        sweep_num = f"Sweep {disp_sweep+1}" if disp_sweep >= 0 else "All Sweeps"
        fdate = f"{fname_parts[-3]}-{fname_parts[-2]}-{fname_parts[-1]} {fname_parts[0]} {sweep_num}"
        fig.suptitle(f"{fdate}")
        ax1.plot(x, y_raw[:, disp_sweep])
        ax2.plot(x, y_processed[:, disp_sweep])

        ax2.set_xlabel("Time (s)")
        ax1.set_ylabel("Raw LFP (V)")
        ax2.set_ylabel("Processed LFP\n(V relative to mean)")

    #plt.xlabel(w['dimension_units'].decode().strip('\x00'))
    #plt.ylabel(w['data_units'].decode().strip('\x00'))

    plt.show()


def plot_scatter_columns(expt_types : list, transients : dict):
    # Plot bar chart of means and CVs by experiment type, before and after PEP introduction
    fig, (ax1) = plt.subplots(1, 1, sharex=True)

    # Check for columns with no entries
    no_entry_keys = []
    max_col_num = len(expt_types)
    col_num = max_col_num
    btbr_max_col = col_num / 2
    for i, k in enumerate(transients.keys()):
        if len(transients[k]) == 0:
            no_entry_keys.append((i, k))
            # Shrink the starting position of the prepep region as it decreases in size
            if i < col_num:
                col_num -= 1
    # Delete all entryless columns
    for i, k in no_entry_keys:
        del transients[k]
        if i < max_col_num:
            del expt_types[i]
        # Keep track of the last BTBR column
        if i < btbr_max_col:
            btbr_max_col -= 1

    # Generate column positions
    col_ids = np.arange(col_num)

    # Construct lists of keys, scatter x positions, and scatter y positions for prepep and postpep
    prepep_keys = list(transients.keys())[:col_num]
    postpep_keys = list(transients.keys())[col_num:]
    prepep_x_list = [col for col, k in enumerate(prepep_keys)
                     for _ in transients[k]]
    prepep_y_list = [v for k in prepep_keys for v in transients[k]]
    postpep_x_list = [col for col, k in enumerate(postpep_keys)
                      for _ in transients[k]]
    postpep_y_list = [v for k in postpep_keys for v in transients[k]]

    # Generate column labels
    col_labels = []
    for t in expt_types:
        label = ""
        label += "BTBR" if "btbr" in t else "WT"
        if "control" in t:
            label += " Control"
        if "cbd" in t:
            label += " CBD"
        if "picro" in t:
            label += " Picrotoxin"
        if "bc" in t:
            label += " BC"
        col_labels.append(label)

    # Plot column scatter plot of transient counts
    ax1.scatter(prepep_x_list, prepep_y_list, color=(0.0, 0.0, 1.0, 0.5),
                marker='o', s=80, label="Pre-PEP")
    ax1.scatter(postpep_x_list, postpep_y_list, color=(1.0, 0.0, 0.0, 0.5),
                marker='o', s=80, label="Post-PEP")
    ax1.set_ylabel("# of Transients", fontsize=20)
    #ax2.set_ylabel("CV of Transients")
    ax1.set_xlabel("Experiment Type", fontsize=20)
    ax1.tick_params(axis='both', which='major', labelsize=15)
    if btbr_max_col > 0:
        ax1.axvline(btbr_max_col - 0.5, c='k')
    ax1.grid(visible=True, axis='y', which='major')

    plt.xticks(col_ids, col_labels)
    plt.legend()
    #ax2.bar(transients.keys(), cv_list.values())
    plt.show()
    #print(f"\nMean coefficient of variation: {np.mean(np.asarray(cv_list))}")