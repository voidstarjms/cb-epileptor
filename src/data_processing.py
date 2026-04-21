import os
import numpy as np
import datetime
import pickle
from brian2 import *
from brian2tools import *
from typing import Dict, Any

_OUTPUT_DATA_FILE = 'output.pkl'


def create_spike_matrix_histo(params_dict: Dict, spike_data: Dict, num_cells: int) -> np.ndarray:
    spike_times = spike_data['t']
    neuron_indices = spike_data['i']

    duration = params_dict['SIM_DURATION'] / second + params_dict['TRANSIENT']
    dt = 0.02  # 20ms per bin

    time_bins = np.arange(params_dict['TRANSIENT'], duration + dt, dt)
    neuron_bins = np.arange(0, num_cells + 1)

    spike_matrix, neuron_edges, time_edges = np.histogram2d(
        neuron_indices,
        spike_times,
        bins=[neuron_bins, time_bins]
    )

    return spike_matrix


def save_data(filepaths: Any, params_dict: Dict, M_N1: Any, M_N2: Any,
              SM_N1: Any, SM_N2: Any, M_S1_1: Any = None, cb_on: bool = True) -> None:
    os.makedirs(filepaths.data_dir, exist_ok=True)
    sim_data = {
        'metadata': {
            'timestamp': datetime.datetime.now().isoformat(),
            'brian2_version': '2.x',
        },
        'params': params_dict,
        'results': {
            't': np.asarray(M_N1.t),
            # POP 1
            'x1': np.asarray(M_N1.x),
            'y1': np.asarray(M_N1.y),
            'z1': np.asarray(M_N1.z),
            'I_syn_inter_1': np.asarray(M_N1.I_syn_inter),
            'I_syn_intra_1': np.asarray(M_N1.I_syn_intra),
            # POP 2
            'x2': np.asarray(M_N2.x),
            'n2': np.asarray(M_N2.n),
            'I_syn_inter_2': np.asarray(M_N2.I_syn_inter),
            # SPIKES
            'spikes_n1': {'t': np.asarray(SM_N1.t), 'i': np.asarray(SM_N1.i)},
            'spikes_n2': {'t': np.asarray(SM_N2.t), 'i': np.asarray(SM_N2.i)},
        }
    }

    if cb_on:
        sim_data['results'].update({
            'syn_wpre': np.asarray(M_S1_1.Wpre),
            'u': np.asarray(M_S1_1.u),
            'Ca': np.asarray(M_S1_1.Ca),
        })

    filepath = os.path.join(filepaths.data_dir, _OUTPUT_DATA_FILE)
    with open(filepath, 'wb') as f:
        pickle.dump(sim_data, f)

    print(f"Simulation data and parameters saved to: {filepath}")


def load_sim_data(filepaths: Any) -> Dict:
    filepath = os.path.join(filepaths.data_dir, _OUTPUT_DATA_FILE)
    with open(filepath, 'rb') as f:
        data = pickle.load(f)
    return data


def dump_spikes_to_file(filename: str, neuron_idx: np.ndarray, spike_times: np.ndarray) -> None:
    """Dump spike times for neuron 0 to a text file for manual correctness checks."""
    mask = np.where(neuron_idx == 0)
    n0_spikes = spike_times[mask]
    np.savetxt(filename, n0_spikes, fmt='%f', delimiter=' ')


def dump_array_to_file(filename: str, arr: np.ndarray) -> None:
    np.savetxt(filename, arr, fmt='%f', delimiter=' ')
