import os
import numpy as np
import config
import params 
import datetime 
import pickle
from brian2 import *
from brian2tools import *

DATA_DIR = config.DATA_DIR
FIGURES_DIR = config.FIGURES_DIR
OUTPUT_DATA_FILE = config.OUTPUT_DATA_FILE

def create_spike_matrix_histo(spike_data, num_cells, transient):
    spike_times = spike_data['t'] 
    neuron_indices = spike_data['i'] 

    duration = params.SIM_DURATION/second
    dt = 0.1  # 100ms per bin

    valid = spike_times > transient
    spike_times = spike_times[valid]
    neuron_indices = neuron_indices[valid]

    time_bins = np.arange(0, duration + dt, dt)
    neuron_bins = np.arange(0, num_cells + 1)

    spike_matrix, neuron_edges, time_edges = np.histogram2d(
        neuron_indices, 
        spike_times,   
        bins=[neuron_bins, time_bins]
    )

    return spike_matrix

def get_params_dict():
    """
    Extracts explicit parameters from the params module.
    Filters out modules and private variables, keeping only uppercase config variables
    that are safe to pickle (scalars, arrays, and Brian2 Quantities).
    """
    params_dict = {}
    for key, val in vars(params).items():
        if key.isupper():
            # Only save specific types. 
            # Brian2 internal structures like DEFAULT_FUNCTIONS cause pickle errors
            if isinstance(val, (int, float, str, bool, np.ndarray, Quantity)):
                params_dict[key] = val
    return params_dict

def save_data(M_N1, M_N2, SM_N1, SM_N2):

    # Save Data, Metadata, and Parameters
    if not os.path.exists(DATA_DIR):
        os.makedirs(DATA_DIR)
    sim_data = {
        'metadata': {
            'timestamp': datetime.datetime.now().isoformat(),
            'brian2_version': '2.x',
            'sim_duration': params.SIM_DURATION
        },
        'params': get_params_dict(),
        'results': {
            't': np.asarray(M_N1.t),
            'x1': np.asarray(M_N1.x),
            'y1': np.asarray(M_N1.y),
            'z1': np.asarray(M_N1.z),
            'I_syn_inter_1': np.asarray(M_N1.I_syn_inter),
            'x2': np.asarray(M_N2.x),
            'n2': np.asarray(M_N2.n),
            'I_syn_inter_2': np.asarray(M_N2.I_syn_inter),
            'spikes_n1': {'t': np.asarray(SM_N1.t), 'i': np.asarray(SM_N1.i)},
            'spikes_n2': {'t': np.asarray(SM_N2.t), 'i': np.asarray(SM_N2.i)}
        }
    }

    # Dump to pickle
    filepath = os.path.join(DATA_DIR, OUTPUT_DATA_FILE)
    with open(filepath, 'wb') as f:
        pickle.dump(sim_data, f)
    
    print(f"Simulation data and parameters saved to: {filepath}")

def load_sim_data():
        # Load pickle
    filepath = os.path.join(DATA_DIR, OUTPUT_DATA_FILE)
    with open(filepath, 'rb') as f:
        data = pickle.load(f)

    return data


def cutoff_transient(data, transient, dt):
    if transient > 0:
        start_idx = int(np.ceil(transient/dt)-1)
        if len(data.shape) == 2:
            data = data[:,start_idx:]
        elif len(data.shape) == 1:
            data = data[start_idx:]
    return data
    

def dump_spikes_to_file(neuron_idx, spike_times):
    mask = np.where(neuron_idx == 0)
    n0_spikes = spike_times[mask]
    
    np.savetxt('spike_times.txt', n0_spikes, fmt='%f', delimiter=' ')

    