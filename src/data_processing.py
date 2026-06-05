"""Save and load sim output, plus a few post-processing helpers."""
import os
import numpy as np
import numpy.typing as npt
import datetime
import pickle
from brian2 import *
from typing import Dict, Any

_OUTPUT_DATA_FILE = 'output.pkl'


def create_spike_matrix_histo(params_dict: Dict, spike_data: Dict, num_cells: int) -> np.ndarray:
    """Bin spikes into a (num_cells, time_bins) count matrix using 20 ms bins.

    Args:
        params_dict (Dict): Needs SIM_DURATION and TRANSIENT.
        spike_data (Dict): Dict with 't' (spike times in seconds) and 'i' (neuron indices).
        num_cells (int): Number of neurons; sets the neuron-axis size.

    Returns:
        np.ndarray: 2D array of shape (num_cells, num_time_bins) holding spike counts.
    """
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


def build_sim_dict(params_dict: Dict, M_N1: Any, M_N2: Any,
                   SM_N1: Any, SM_N2: Any, M_S1_1: Any = None,
                   cb_on: bool = True, M_param: Any = None) -> Dict:
    """Build the sim_data dict from monitors (no I/O).

    Args:
        params_dict (Dict): Parameter dict from param_loader; embedded in output.
        M_N1 (Any): StateMonitor on the HR population.
        M_N2 (Any): StateMonitor on the ML population.
        SM_N1 (Any): SpikeMonitor on the HR population.
        SM_N2 (Any): SpikeMonitor on the ML population.
        M_S1_1 (Any): Optional StateMonitor on HR->HR synapses (Wpre, u, Ca).
        cb_on (bool): If True, include plasticity traces (Wpre, u, Ca) in the output.
        M_param (Any): Optional StateMonitor recording the x0/ce schedule traces.

    Returns:
        Dict: {'metadata', 'params', 'results'} ready to pickle or pass downstream.
    """
    sim_data = {
        'metadata': {
            'timestamp': datetime.datetime.now().isoformat(),
            'brian2_version': '2.x',
        },
        'params': params_dict,
        'results': {
            't': np.asarray(M_N1.t),
            'x1': np.asarray(M_N1.x),
            'x2': np.asarray(M_N2.x),
            'I_syn_inter_2': np.asarray(M_N2.I_syn_inter),
            'spikes_n1': {'t': np.asarray(SM_N1.t), 'i': np.asarray(SM_N1.i)},
            'spikes_n2': {'t': np.asarray(SM_N2.t), 'i': np.asarray(SM_N2.i)},
        }
    }

    if cb_on and M_S1_1 is not None:
        sim_data['results'].update({
            'syn_wpre': np.asarray(M_S1_1.Wpre),
            'u':        np.asarray(M_S1_1.u),
            'Ca':       np.asarray(M_S1_1.Ca),
        })

    if M_param is not None:
        sim_data['results'].update({
            'x0_t':  np.asarray(M_param.x0_t[0]),
            'ce_t':  np.asarray(M_param.ce_t[0]),
        })

    return sim_data


def save_sim_dict(filepaths: Any, sim_data: Dict) -> None:
    """Pickle a sim_data dict to filepaths.data_dir/output.pkl.

    Args:
        filepaths (Any): FilePaths with data_dir.
        sim_data (Dict): Dict from build_sim_dict.
    """
    os.makedirs(filepaths.data_dir, exist_ok=True)
    filepath = os.path.join(filepaths.data_dir, _OUTPUT_DATA_FILE)
    with open(filepath, 'wb') as f:
        pickle.dump(sim_data, f)

    print(f"Simulation data and parameters saved to: {filepath}")


def load_sim_data(filepaths: Any) -> Dict:
    """Return the dict written by save_sim_dict.

    Args:
        filepaths (Any): FilePaths with data_dir.

    Returns:
        Dict: Dict with keys 'metadata', 'params', 'results'.
    """
    filepath = os.path.join(filepaths.data_dir, _OUTPUT_DATA_FILE)
    with open(filepath, 'rb') as f:
        data = pickle.load(f)
    return data

def delete_sim_data(filepaths: Any) -> None:
    """Delete the dict written by save_sim_dict.

    Args:
        filepaths (Any): FilePaths dataclass with data_dir.
    """
    filepath = os.path.join(filepaths.data_dir, _OUTPUT_DATA_FILE)
    os.remove(filepath)

def dump_spikes_to_file(filename: str, neuron_idx: np.ndarray, spike_times: np.ndarray) -> None:
    """Dump spike times for neuron 0 to a text file for manual correctness checks.

    Args:
        filename (str): Output text file path.
        neuron_idx (np.ndarray): Per-spike neuron index (e.g. SpikeMonitor.i).
        spike_times (np.ndarray): Per-spike times (e.g. SpikeMonitor.t).
    """
    mask = np.where(neuron_idx == 0)
    n0_spikes = spike_times[mask]
    np.savetxt(filename, n0_spikes, fmt='%f', delimiter=' ')


def dump_array_to_file(filename: str, arr: np.ndarray) -> None:
    """Write arr as whitespace-separated floats.

    Args:
        filename (str): Output text file path.
        arr (np.ndarray): Array to write.
    """
    np.savetxt(filename, arr, fmt='%f', delimiter=' ')



def save_raw_spikes(filepaths: Any, spike_data: Dict, filename: str = "sparse_spike_data.pkl") -> None:
    """Test function to benchmark data save sizes. This function is less memory efficient and not used"""
    os.makedirs(filepaths.data_dir, exist_ok=True)
    print(f"+++TYPE t: {spike_data['t'].dtype}, i: {spike_data['i'].dtype}+++")
    save_spike_data = {
        "t": spike_data["t"],
        "i": spike_data["i"]
    }

    filepath = os.path.join(filepaths.data_dir, filename)
    with open(filepath, 'wb') as f:
        pickle.dump(save_spike_data, f)

    print(f"Simulation data and parameters saved to: {filepath}")

def save_spike_histo(filepaths: Any, spike_matrix: npt.NDArray[np.float64], filename: str = "spike_matrix_data.pkl") -> None:
    """Test function to benchmark data save sizes. This function is the better method. Still not used, but run_single_sim does a similar thing as this"""
    os.makedirs(filepaths.data_dir, exist_ok=True)
    spike_matrix = spike_matrix.astype(np.uint8)
    print(f"+++TYPE: {spike_matrix.dtype}+++")
    save_spike_data = {
        "spikes": spike_matrix
    }
    
    filepath = os.path.join(filepaths.data_dir, filename)
    with open(filepath, 'wb') as f:
        pickle.dump(save_spike_data, f)

    print(f"Simulation data and parameters saved to: {filepath}")