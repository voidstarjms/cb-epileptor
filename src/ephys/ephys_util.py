import numpy as np
from scipy.stats import wilcoxon

def separate_prepost_columns(expt_types : list, transients : dict):
    """Separate transient count lists by pre vs. post.
    Delete experiment types with empty transient count lists.
    Compute split point for the two arms of the experiments.

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
            # Keep track of the first column of
            # experiment arm 2 (e.g. BTBR)
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

    return (expt_types, col_num, split_point)

def classwise_wilcoxon(transients, expt_types):
    """Compute Wilcoxon signed-rank test for each experiment class (pre-PEP against post-PEP)
    
    Args:
        transients (dict[list]): Dictionary of lists containing transient counts for each
            experiment of each experiment type.
        expt_types (list[str]): List of key strings for each experiment type.
    
    Returns:
        tuple (stat_list, pvalue_list):
            stat_list: statistic values for each experiment type.
            pvalue_list: p-values for each experiment type.
    """
    stat_list = []
    pvalue_list = []
    prepep_keys = [ent+"_pre" for ent in expt_types]
    postpep_keys = [ent+"_post" for ent in expt_types]
    for pre_k, post_k in zip(prepep_keys, postpep_keys):
        x = np.array(transients[pre_k])
        y = np.array(transients[post_k])
        wilcoxon_results = wilcoxon(y - x)
        stat_list.append(wilcoxon_results.statistic)
        pvalue_list.append(wilcoxon_results.pvalue)

    return stat_list, pvalue_list