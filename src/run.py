import argparse
import os
import shutil
import numpy as np
from brian2 import *
from dataclasses import dataclass
import data_processing
import synch as syn
import plotting.population_plots as pop_plotter
import plotting.plasticity_plots as plast_plotter
import plotting.analysis_plots as analysis_plotter
from model import run_sim
from typing import Dict
import yaml


@dataclass
class FilePaths:
    """Groups output paths; instantiated from CLI args in main() and injected throughout."""
    data_dir: str
    figures_dir: str


def plot_output(filepaths: FilePaths, params_dict: Dict, cb_on: bool = True) -> None:
    os.makedirs(filepaths.figures_dir, exist_ok=True)
    data = data_processing.load_sim_data(filepaths)
    res = data['results']
    t = res['t']
    x1 = res['x1']
    x2 = res['x2']

    if cb_on:
        wpre = res['syn_wpre'][0]
        u = res['u'][0]
        Ca = res['Ca'][0]
        plast_plotter.plot_wpre(filepaths, params_dict, t, x1, wpre, u, Ca)

    num_cells = params_dict['NUM_CELLS']
    spike_matrix_1 = data_processing.create_spike_matrix_histo(params_dict, res['spikes_n1'], num_cells)
    spike_matrix_2 = data_processing.create_spike_matrix_histo(params_dict, res['spikes_n2'], num_cells)
    pop_plotter.standard_plot(filepaths, params_dict, t, x1, x2, spike_matrix_1, spike_matrix_2,
                              num_cells, params_dict['SIM_DURATION'] / second,
                              g_inter_vals=params_dict['G_INTER_VALS'], g_intra_vals=params_dict['G_INTRA_VALS'],
                              coupling_vals=params_dict['COUPLING_VALS'], x_naught_vals=params_dict['X_NAUGHT_VALS'])


def plot_output_full(filepaths: FilePaths, params_dict: Dict, cb_on: bool = True) -> None:
    os.makedirs(filepaths.figures_dir, exist_ok=True)
    data = data_processing.load_sim_data(filepaths)
    res = data['results']
    t = res['t']
    x1 = res['x1']
    y1 = res['y1']
    z1 = res['z1']
    I_syn_inter = res['I_syn_inter_1']
    I_syn_intra = res['I_syn_intra_1']
    x2 = res['x2']
    n = res['n2']

    if cb_on:
        wpre = res['syn_wpre'][0]
        u = res['u'][0]
        Ca = res['Ca'][0]
        plast_plotter.plot_wpre(filepaths, params_dict, t, x1, wpre, u, Ca)

    num_cells = params_dict['NUM_CELLS']
    spike_matrix_1 = data_processing.create_spike_matrix_histo(params_dict, res['spikes_n1'], num_cells)
    spike_matrix_2 = data_processing.create_spike_matrix_histo(params_dict, res['spikes_n2'], num_cells)

    pop_plotter.plot_hr_single(filepaths, t, x1, y1, z1, I_syn_inter)
    pop_plotter.plot_ml_single(filepaths, t, x2, n)
    pop_plotter.standard_plot(filepaths, params_dict, t, x1, x2, spike_matrix_1, spike_matrix_2,
                              num_cells, params_dict['SIM_DURATION'] / second,
                              g_inter_vals=params_dict['G_INTER_VALS'], g_intra_vals=params_dict['G_INTRA_VALS'],
                              coupling_vals=params_dict['COUPLING_VALS'], x_naught_vals=params_dict['X_NAUGHT_VALS'])

    print("I_syn_inter max (raw amps):", np.max((I_syn_inter[0] / amp)))
    print("I_syn_intra max (raw amps):", np.max((I_syn_intra[0] / amp)))
    print("I_syn_inter min (raw amps):", np.min((I_syn_inter[0] / amp)))
    print("I_syn_intra min (raw amps):", np.min((I_syn_intra[0] / amp)))


def analyze_populations(filepaths: FilePaths, params_dict: Dict, data: Dict) -> None:
    res = data['results']
    x1 = res['x1']
    x2 = res['x2']
    pop1_spike_times, pop1_neuron_idx = res['spikes_n1']['t'], res['spikes_n1']['i']
    pop2_spike_times, pop2_neuron_idx = res['spikes_n2']['t'], res['spikes_n2']['i']

    data_processing.dump_spikes_to_file(
        os.path.join(filepaths.data_dir, 'spikes_n1.txt'),
        np.asarray(pop1_neuron_idx),
        np.asarray(pop1_spike_times),
    )

    print("============HINDMARSH ROSE STATS============")
    chi, autocorr, lag = syn.autocorrelate(x1)
    print(f'synchrony measure: {chi}\nautocorrelation: {autocorr}')
    z, r, psi = syn.KOP(pop1_neuron_idx, pop1_spike_times, params_dict['SIM_DURATION'] / second)
    print(f'r: {np.mean(r)}')
    analysis_plotter.plot_autocorr(filepaths, autocorr, lag)
    phase_matrix = syn.compute_phase(pop1_neuron_idx, pop1_spike_times, params_dict['SIM_DURATION'] / second)
    analysis_plotter.plot_kop(filepaths, phase_matrix)

    print("\n============MORRIS LECAR STATS============")
    chi, autocorr, lag = syn.autocorrelate(x2)
    print(f'synchrony measure: {chi}\nautocorrelation: {autocorr}')
    z, r, psi = syn.KOP(pop2_neuron_idx, pop2_spike_times, params_dict['SIM_DURATION'] / second)
    print(f'r: {np.mean(r)}')


def _save_params(params_file: str, run_dir: str) -> None:
    shutil.copy2(params_file, os.path.join(run_dir, os.path.basename(params_file)))


REQUIRED_PARAMS = {
    'SIM_DURATION', 'NUM_CELLS', 'TAU_CLOCK', 'DT_SCALING', 'TRANSIENT',
    'ISOLATE', 'W_MAX', 'I_SCALE', 'NOISE_INIT_OFFSET',
    'HR_A', 'HR_B', 'HR_C', 'HR_D', 'HR_S', 'HR_I_APP',
    'HR_X_NAUGHT', 'HR_R', 'HR_SIGMA', 'HR_THRESHOLD', 'HR_REFRACTORY_CONDITION',
    'ML_CM', 'ML_I_APP', 'ML_GL', 'ML_E_L', 'ML_GK', 'ML_E_K',
    'ML_GCA', 'ML_E_CA', 'ML_V1', 'ML_V2', 'ML_V3', 'ML_V4',
    'ML_PHI', 'ML_SIGMA', 'ML_Z_BAR_SCALE', 'ML_Z_BAR_OFFSET',
    'ML_THRESHOLD', 'ML_REFRACTORY_CONDITION',
    'SYN_TMAX', 'SYN_VT', 'SYN_KP',
    'SYN_ALPHA_EXC', 'SYN_BETA_EXC', 'SYN_E_EXC',
    'SYN_ALPHA_INH', 'SYN_BETA_INH', 'SYN_E_INH',
    'THETA_LTD_START', 'THETA_LTD_END', 'THETA_LTP_START',
    'A_LTD', 'A_LTP', 'TAU_WPRE', 'TAU_CA',
    'CA_SIGMOID_SHIFT', 'CA_SIGMOID_SLOPE',
    'X_NAUGHT_VALS', 'COUPLING_VALS', 'G_INTER_VALS', 'G_INTRA_VALS',
}

# Maps YAML param keys to their Brian2 units. Keys absent from this dict are unitless.
_PARAM_UNITS: Dict = {
    'SIM_DURATION':  second,
    'TAU_CLOCK':     msecond,
    'TAU_WPRE':      second,
    'TAU_CA':        msecond,
    'I_SCALE':       uamp,
    'HR_R':          1 / msecond,
    'ML_CM':         ufarad,
    'ML_I_APP':      uamp,
    'ML_V1':  mvolt, 'ML_V2': mvolt, 'ML_V3': mvolt, 'ML_V4': mvolt,
    'ML_PHI':        1 / msecond,
    'ML_E_CA': mvolt, 'ML_E_K': mvolt, 'ML_E_L': mvolt,
    'ML_GL':  msiemens, 'ML_GCA': msiemens, 'ML_GK': msiemens,
    'ML_SIGMA':      uA,
    'SYN_VT':  mV,   'SYN_KP': mV,
    'SYN_TMAX':      mmolar,
    'SYN_ALPHA_EXC': 1 / (mmolar * msecond),
    'SYN_BETA_EXC':  1 / msecond,
    'SYN_E_EXC':     mV,
    'SYN_ALPHA_INH': 1 / (mmolar * msecond),
    'SYN_BETA_INH':  1 / msecond,
    'SYN_E_INH':     mV,
    'G_INTER_VALS':  uS,
    'G_INTRA_VALS':  uS,
}


def load_params(params_file: str) -> Dict:
    """Load and validate a YAML params file."""
    ext = os.path.splitext(params_file)[1].lower()
    if ext not in ('.yaml', '.yml'):
        raise ValueError(f"Unsupported params file type: {ext}. Expected .yaml or .yml")
    params_dict = _load_params_yaml(params_file)
    missing = REQUIRED_PARAMS - params_dict.keys()
    if missing:
        raise ValueError(f"Missing required params in {params_file}: {sorted(missing)}")
    return params_dict


def _load_params_yaml(params_file: str) -> Dict:
    with open(params_file, 'r') as f:
        raw = yaml.safe_load(f)
    return {
        k: (v * _PARAM_UNITS[k] if k in _PARAM_UNITS else v)
        for k, v in raw.items()
    }


def main() -> None:
    DEFAULT_OUT_DIR = 'output/'
    DEFAULT_PARAMS = 'parameters/params1.yaml'

    parser = argparse.ArgumentParser(description="Run and/or plot the simulation.")
    parser.add_argument('-m', '--mode', type=str, default='rp',
                        help="Run mode: 'r' run, 'p' plot, 'a' analyze, 't' test.")
    parser.add_argument('--cb', action='store_true', default=True,
                        help="Enable CB synapses (default: enabled)")
    parser.add_argument('--no-cb', dest='cb', action='store_false',
                        help="Disable CB synapses")
    parser.add_argument('--params', type=str, default=DEFAULT_PARAMS,
                        help="Parameter set to use (default: 'default')")
    parser.add_argument('--out-dir', type=str, default=DEFAULT_OUT_DIR,
                        help="Output directory (default: 'output/')")
    args = parser.parse_args()
    run_mode = args.mode
    cb_on = args.cb
    params = args.params
    params_dict = load_params(params)

    out_dir = args.out_dir
    filepaths = FilePaths(
        data_dir=os.path.join(out_dir, 'data'),
        figures_dir=os.path.join(out_dir, 'figures'),
    )

    if 'r' in run_mode:
        os.makedirs(filepaths.data_dir)
        os.makedirs(filepaths.figures_dir)
        _save_params(params, out_dir)
        print("Running simulation...")
        run_sim(filepaths, params_dict, cb_on)
        print(f"Simulation complete. Results saved to {out_dir}")

    if 'p' in run_mode:
        print("Generating plots...")
        if 'f' in run_mode:
            plot_output_full(filepaths, params_dict, cb_on)
        else:
            plot_output(filepaths, params_dict, cb_on)
        print(f"Plots saved to {filepaths.figures_dir}.")

    if 'a' in run_mode:
        data = data_processing.load_sim_data(filepaths)
        analyze_populations(filepaths, data['params'], data)


if __name__ == "__main__":
    main()
