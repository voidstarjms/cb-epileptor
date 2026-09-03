import numpy as np
from scipy.stats import wilcoxon
import sys

EXPT_TYPE_LIST = None
PREPEP_KEYS = None
POSTPEP_KEYS = None

def separate_prepost_columns(expt_types : list, transients : dict):
    """Separate transient count lists by pre vs. post.
    Delete experiment types with empty transient count lists.
    Compute split point for the two arms of the experiments.

    Args:
        expt_types (list[str]): List of experiment type names, keys to the transients dict.
        transients (dict[list]): Dictionary of lists containing transient counts for each
            segment of each experiment of each experiment type.
        
    Returns:
        tuple:
            expt_types (list[str]): The original expt_types list with entries which are keys
                to empty transient lists removed.
            split_point (int): The column at which the second arm of the experiment type set
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

    return (expt_types, split_point)

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

def get_expt_types(type_file_path=None):
    global EXPT_TYPE_LIST # I hate this

    # Initialize the singleton
    if EXPT_TYPE_LIST == None:
        if type_file_path == None:
            print("get_expt_types: singleton is not initialized and no file was specified.")
            sys.exit(1)

        EXPT_TYPE_LIST = []
        control = False
        cbd = False
        picro = False
        cbd_picro = False
        cbd_bc = False
        cbd_amounts = None
        cbd_list = []

        # Read CSV to get experiment types and CBD concentrations
        with open(type_file_path, "r") as f:
            expt_types = f.readline().split(",")
            if len(expt_types) == 0:
                print("_init_transient_dict: expt_types.csv is empty or starts with empty line.")
                sys.exit(1)
            for i in range(len(expt_types)):
                expt_types[i] = expt_types[i].strip()
            
            for ent in expt_types:
                match ent:
                    case "control":
                        control = True
                    case "cbd":
                        cbd = True
                    case "picro":
                        picro = True
                    case "cbd_picro":
                        cbd_picro = True
                    case "cbd_bc":
                        cbd_bc = True

            cbd_amounts = f.readline().split(",")
            if cbd_amounts != None:
                for i in range(len(cbd_amounts)):
                    cbd_amounts[i] = cbd_amounts[i].strip()

        # Generate prefixes for all CBD concentrations
        if cbd_amounts is None:
            cbd_list.append("cbd")
        else:
            for i in cbd_amounts:
                cbd_list.append("cbd_"+i)

        # Control entries
        if control == True:
            EXPT_TYPE_LIST.append("control")
        
        # CBD entries
        if cbd == True:
            for i in cbd_list:
                EXPT_TYPE_LIST.append(i)

        # Picro
        if picro == True:
            EXPT_TYPE_LIST.append("picro")

        # CBD picro
        if cbd_picro == True:
            for i in cbd_list:
                EXPT_TYPE_LIST.append(i+"_picro")

        # CBD BC
        if cbd_bc == True:
            for i in cbd_list:
                EXPT_TYPE_LIST.append(i+"_bc")

        # WT vs. BTBR experiments
        EXPT_TYPE_LIST = ["wt_"+t for t in EXPT_TYPE_LIST] +\
            ["btbr_"+t for t in EXPT_TYPE_LIST]
    
    return EXPT_TYPE_LIST

def get_expt_keys(type_file_path=None):
    global PREPEP_KEYS
    global POSTPEP_KEYS

    # Initialize the singletons
    if PREPEP_KEYS == None:
        get_expt_types(type_file_path)
        PREPEP_KEYS = [ent+"_pre" for ent in EXPT_TYPE_LIST]
        POSTPEP_KEYS = [ent+"_post" for ent in EXPT_TYPE_LIST]
    
    return PREPEP_KEYS, POSTPEP_KEYS

def split_prepost_transients(transients, type_file_path=None):
    global PREPEP_KEYS
    global POSTPEP_KEYS

    get_expt_keys(type_file_path)
    pre = [e for k in PREPEP_KEYS for e in transients[k]]
    post = [e for k in POSTPEP_KEYS for e in transients[k]]
    return pre, post