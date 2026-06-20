import argparse
import os
import shutil
import numpy as np
from brian2 import *
from dataclasses import dataclass
import data_processing
import synch_metrics as syn
import plotting.population_plots as pop_plotter
import plotting.plasticity_plots as plast_plotter
import plotting.analysis_plots as analysis_plotter
from model import run_sim
from param_loader import load_params
from config_structs import ZoomConfig
from typing import Dict


@dataclass
class FilePaths:
    """Groups output paths; instantiated from CLI args in main() and injected throughout.

    Attributes:
        data_dir (str): Directory for sim output (output.pkl, spike dumps).
        figures_dir (str): Directory for generated plots.
    """
    data_dir: str
    figures_dir: str

def plot_synapse_dynamics(filepaths: FilePaths, params_dict: Dict, res: Dict):
    """Plot mean of synapse parameters Wpre, u, and Ca for each population.
    
    Args:
        filepaths (FilePaths): Reads from data_dir, writes to figures_dir.
        params_dict (Dict): Parameter dict loaded from YAML.
        res: Results from the simulation dictionary.
    """
    t = res['t']
    x1 = res['x1']
    x2 = res['x2']

    # E->E synapse stats
    wpre = res['S1_1_wpre']
    u = res['S1_1_u']
    Ca = res['S1_1_Ca']
    plast_plotter.plot_plasticity(filepaths, "S1_1_synapse_stats.png", params_dict, "E-to-E Synapse Stats", t, x1, wpre, u, Ca)
    # E->I synapse stats
    wpre = res['S1_2_wpre']
    u = res['S1_2_u']
    Ca = res['S1_2_Ca']
    plast_plotter.plot_plasticity(filepaths, "S1_2_synapse_stats.png", params_dict, "E-to-I Synapse Stats", t, x1, wpre, u, Ca)
    # I->I synapse stats
    wpre = res['S2_2_wpre']
    u = res['S2_2_u']
    Ca = res['S2_2_Ca']
    plast_plotter.plot_plasticity(filepaths, "S2_2_synapse_stats.png", params_dict, "I-to-I Synapse Stats", t, x2, wpre, u, Ca)
    # I->E synapse stats
    wpre = res['S2_1_wpre']
    u = res['S2_1_u']
    Ca = res['S2_1_Ca']
    plast_plotter.plot_plasticity(filepaths, "S2_1_synapse_stats.png", params_dict, "I-to-E Synapse Stats", t, x2, wpre, u, Ca)


def plot_output(filepaths: FilePaths, params_dict: Dict, cb_on: bool = True, results_dict: Dict = None) -> None:
    """Load sim data and make the default plots (LFP + both rasters).

    Args:
        filepaths (FilePaths): Reads from data_dir, writes to figures_dir.
        params_dict (Dict): Parameter dict loaded from YAML.
        cb_on (bool): If True, also plot the plasticity (Wpre, u, Ca) traces.
        results_dict (Dict): Run result dictionary, used if not None.
    """
    
    os.makedirs(filepaths.figures_dir, exist_ok=True)
    res = results_dict['results']
    t = res['t']
    x1 = res['x1']
    x2 = res['x2']

    if cb_on:
        plot_synapse_dynamics(filepaths, params_dict, res)
    
    num_cells = params_dict['NUM_CELLS']
    spike_matrix_1 = data_processing.create_spike_matrix_histo(params_dict, res['spikes_n1'], num_cells)
    spike_matrix_2 = data_processing.create_spike_matrix_histo(params_dict, res['spikes_n2'], num_cells)
    pop_plotter.standard_plot(filepaths, params_dict, t, x1, x2, spike_matrix_1, spike_matrix_2,
                              num_cells, params_dict['SIM_DURATION'] / second,
                              g_inter_vals=params_dict['G_INTER_VALS'], g_intra_vals=params_dict['G_INTRA_VALS'],
                              x0_t=res.get('x0_t'), ce_t=res.get('ce_t'))

# TODO This is presently deprecated because we're not saving y, z, etc. Figure out what to do with it
def plot_output_full(filepaths: FilePaths, params_dict: Dict, cb_on: bool = True) -> None:
    """Same as plot_output but also plots EPOP (x, y, z, I_syn_inter) and IPOP (x, n) traces.

    Args:
        filepaths (FilePaths): Reads from data_dir, writes to figures_dir.
        params_dict (Dict): Parameter dict loaded from YAML.
        cb_on (bool): If True, also plot the plasticity (Wpre, u, Ca) traces.
    """
    os.makedirs(filepaths.figures_dir, exist_ok=True)
    data = data_processing.load_sim_data(filepaths)
    res = data['results']
    t = res['t']
    x1 = res['x1']
    # y1 = res['y1']
    # z1 = res['z1']
    # I_syn_inter = res['I_syn_inter_1']
    # I_syn_intra = res['I_syn_intra_1']
    x2 = res['x2']
    # n = res['n2']

    if cb_on:
        plot_synapse_dynamics(filepaths, params_dict, res)

    num_cells = params_dict['NUM_CELLS']
    spike_matrix_1 = data_processing.create_spike_matrix_histo(params_dict, res['spikes_n1'], num_cells)
    spike_matrix_2 = data_processing.create_spike_matrix_histo(params_dict, res['spikes_n2'], num_cells)

    zoom = ZoomConfig(start=60, end=60.2)
    pop_plotter.plot_hr_single(filepaths, t, x1, y1, z1, I_syn_inter, zoom=zoom)
    pop_plotter.plot_ml_single(filepaths, t, x2, n, zoom=zoom)
    pop_plotter.standard_plot(filepaths, params_dict, t, x1, x2, spike_matrix_1, spike_matrix_2,
                              num_cells, params_dict['SIM_DURATION'] / second, zoom=zoom)
                            #   g_inter_vals=params_dict['G_INTER_VALS'], g_intra_vals=params_dict['G_INTRA_VALS'],
                            #   x0_t=res.get('x0_t'), ce_t=res.get('ce_t'))

    # print("I_syn_inter max (raw amps):", np.max((I_syn_inter[0] / amp)))
    # print("I_syn_intra max (raw amps):", np.max((I_syn_intra[0] / amp)))
    # print("I_syn_inter min (raw amps):", np.min((I_syn_inter[0] / amp)))
    # print("I_syn_intra min (raw amps):", np.min((I_syn_intra[0] / amp)))


def analyze_populations(filepaths: FilePaths, params_dict: Dict, data: Dict) -> None:
    """Print chi and mean KOP r for both pops. Writes autocorr/KOP plots and a spike dump.

    Args:
        filepaths (FilePaths): Spike dump goes to data_dir, plots to figures_dir.
        params_dict (Dict): Parameter dict loaded from YAML.
        data (Dict): Dict returned by data_processing.load_sim_data.
    """
    res = data['results']
    x1 = res['x1']
    x2 = res['x2']
    pop1_spike_times, pop1_neuron_idx = res['spikes_n1']['t'], res['spikes_n1']['i']
    pop2_spike_times, pop2_neuron_idx = res['spikes_n2']['t'], res['spikes_n2']['i']

    # data_processing.dump_spikes_to_file(
    #     os.path.join(filepaths.data_dir, 'spikes_n1.txt'),
    #     np.asarray(pop1_neuron_idx),
    #     np.asarray(pop1_spike_times),
    # )

    # print("===========spike train formats===============")
    # print(pop1_neuron_idx, f"neuron array len is {len(pop1_neuron_idx)}")
    # print(pop1_spike_times, f"spike array len is {len(pop1_spike_times)}")

    num_cells = params_dict['NUM_CELLS']
    data_processing.save_raw_spikes(filepaths, res['spikes_n1'], filename="sparse_spikes_n1.pkl")
    data_processing.save_raw_spikes(filepaths, res['spikes_n2'], filename="sparse_spikes_n2.pkl")
    spike_matrix_1 = data_processing.create_spike_matrix_histo(params_dict, res['spikes_n1'], num_cells)
    spike_matrix_2 = data_processing.create_spike_matrix_histo(params_dict, res['spikes_n2'], num_cells)
    data_processing.save_spike_histo(filepaths, spike_matrix_1, filename="spike_histo_n1.pkl")
    data_processing.save_spike_histo(filepaths, spike_matrix_2, filename="spike_histo_n2.pkl")

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

    analysis_plotter.plot_power_spec(filepaths, params_dict, x1, x2)
    analysis_plotter.plot_spectrogram(filepaths, params_dict, x1, x2,
                                      t=res['t'],
                                      x0_t=res.get('x0_t'),
                                      ce_t=res.get('ce_t'))


def report_metrics(params_dict: Dict, data: Dict) -> None:
    """Print a comparison table of chi, KOP r, mean firing rate, median frequency, and ISI CV.

    Args:
        params_dict: Parameter dict loaded from YAML.
        data: Dict returned by data_processing.load_sim_data.
    """
    res = data['results']
    duration = params_dict['SIM_DURATION'] / second
    num_cells = params_dict['NUM_CELLS']
    fs = float(1 / (params_dict['TAU_CLOCK'] / params_dict['DT_SCALING']) / Hz)

    chi1, _, _ = syn.autocorrelate(res['x1'])
    chi2, _, _ = syn.autocorrelate(res['x2'])

    _, r1, _ = syn.KOP(res['spikes_n1']['i'], res['spikes_n1']['t'], duration)
    _, r2, _ = syn.KOP(res['spikes_n2']['i'], res['spikes_n2']['t'], duration)
    r1, r2 = float(np.mean(r1)), float(np.mean(r2))

    fr1 = syn.mean_firing_rate(res['spikes_n1']['i'], res['spikes_n1']['t'], num_cells, duration)
    fr2 = syn.mean_firing_rate(res['spikes_n2']['i'], res['spikes_n2']['t'], num_cells, duration)

    df1 = syn.median_frequency(res['x1'], fs)
    df2 = syn.median_frequency(res['x2'], fs)

    cv1 = syn.isi_cv(res['spikes_n1']['i'], res['spikes_n1']['t'])
    cv2 = syn.isi_cv(res['spikes_n2']['i'], res['spikes_n2']['t'])

    W = 28
    print("\n============ METRIC REPORT ============")
    print(f"\n{'Metric':<{W}} {'EPOP':>10} {'IPOP':>10}")
    print("-" * (W + 22))
    print(f"{'Chi (autocorr)':<{W}} {chi1:>10.4f} {chi2:>10.4f}")
    print(f"{'KOP r':<{W}} {r1:>10.4f} {r2:>10.4f}")
    print(f"{'Mean firing rate (Hz)':<{W}} {fr1:>10.4f} {fr2:>10.4f}")
    print(f"{'Median frequency (Hz)':<{W}} {df1:>10.4f} {df2:>10.4f}")
    print(f"{'ISI CV':<{W}} {cv1:>10.4f} {cv2:>10.4f}")


def _save_params(params_file: str, run_dir: str) -> None:
    """Copy the params YAML into run_dir so each run is reproducible from its output.

    Args:
        params_file (str): Source YAML path.
        run_dir (str): Destination directory.
    """
    shutil.copy2(params_file, os.path.join(run_dir, os.path.basename(params_file)))


REQUIRED_PARAMS = {
    'SIM_DURATION', 'NUM_CELLS', 'TAU_CLOCK', 'DT_SCALING', 'TRANSIENT',
    'ISOLATE', 'W_MAX', 'I_SCALE', 'NOISE_INIT_OFFSET',
    'EPOP_A', 'EPOP_B', 'EPOP_C', 'EPOP_D', 'EPOP_S', 'EPOP_I_APP',
    'EPOP_BASE_EXCITE', 'EPOP_R', 'EPOP_SIGMA', 'EPOP_THRESHOLD', 'EPOP_REFRACTORY_CONDITION',
    'IPOP_CM', 'IPOP_I_APP', 'IPOP_GL', 'IPOP_E_L', 'IPOP_GK', 'IPOP_E_K',
    'IPOP_GCA', 'IPOP_E_CA', 'IPOP_H_CA', 'IPOP_LAMBDA_CA', 'IPOP_H_K', 'IPOP_LAMBDA_K',
    'IPOP_PHI', 'IPOP_SIGMA', 'IPOP_Z_BAR_SCALE', 'IPOP_Z_BAR_OFFSET',
    'IPOP_THRESHOLD', 'IPOP_REFRACTORY_CONDITION',
    'SYN_TMAX', 'SYN_VT', 'SYN_KP', 'PCT_CONNECT',
    'SYN_ALPHA_EXC', 'SYN_BETA_EXC', 'SYN_E_EXC',
    'SYN_ALPHA_INH', 'SYN_BETA_INH', 'SYN_E_INH',
    'THETA_LTD_START', 'THETA_LTD_END', 'THETA_LTP_START',
    'A_LTD', 'A_LTP', 'TAU_WPRE', 'TAU_CA',
    'CA_SIGMOID_SHIFT', 'CA_SIGMOID_SLOPE',
    'BASE_EXCITE_VALS', 'BASE_EXCITE_DT',
    'COUPLING_VALS', 'COUPLING_DT',
    'G_INTER_VALS', 'G_INTER_DT',
    'G_INTRA_VALS', 'G_INTRA_DT',
}


def main() -> None:
    """Parse CLI args and run the phases selected by --mode.

    Mode flags ('r', 'p', 'a', 'f', 's') are independent and stackable: e.g.
    'rp' runs then plots, 'rpf' runs then makes full plots, 'a' analyzes
    output.
    """
    DEFAULT_OUT_DIR = 'output/run3'
    DEFAULT_PARAMS = 'parameters/params.yaml'   # run.py runs from src/; YAML at repo root

    parser = argparse.ArgumentParser(description="Run and/or plot the simulation.")
    parser.add_argument('-m', '--mode', type=str, default='rp',
                        help="Run mode: 'r' run, 'p' plot, 'a' analyze, 't' test.")
    parser.add_argument('--cb', action='store_true', default=False,
                        help="Enable CB synapses (default: disabled)")
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
    missing = REQUIRED_PARAMS - params_dict.keys()
    if missing:
        raise ValueError(f"Missing required params in {params}: {sorted(missing)}")

    out_dir = args.out_dir
    filepaths = FilePaths(
        data_dir=os.path.join(out_dir, 'data'),
        figures_dir=os.path.join(out_dir, 'figures'),
    )
    
    data = None
    save = 's' in run_mode
    if 'r' in run_mode:
        # Run simulation, save parameters always and results if specified
        os.makedirs(filepaths.data_dir, exist_ok=True)
        os.makedirs(filepaths.figures_dir, exist_ok=True)
        _save_params(params, out_dir)
        print("Running simulation...")
        data = run_sim(filepaths, params_dict, cb_on, save)
        if save:
            print(f"Simulation complete. Results saved to {out_dir}")
    else:
        # Load data from filepaths if simulation wasn't run
        print(f"Loading data from {filepaths.data_dir}")
        data = data_processing.load_sim_data(filepaths)

    # Plot model results
    if 'p' in run_mode:
        print("Generating plots...")
        if 'f' in run_mode:
            print("WARNING: RUN FLAG f IS CURRENTLY DISABLED. SKIPPING PLOTS.")
            #plot_output_full(filepaths, params_dict, cb_on)
        else:
            plot_output(filepaths, params_dict, cb_on, data)
        print(f"Plots saved to {filepaths.figures_dir}.")

    # Analyze model synchrony
    if 'a' in run_mode:
        analyze_populations(filepaths, data['params'], data)

    # Print population behavior metrics
    if 't' in run_mode:
        report_metrics(data['params'], data)


if __name__ == "__main__":
    main()
