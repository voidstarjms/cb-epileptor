import argparse
import os
import numpy as np
from brian2 import second, amp
import params
import config
import data_processing
import synch as syn
import plotting.population_plots as pph
import plotting.plasticity_plots as plph
import plotting.analysis_plots as aph
from simulation.core import run_sim

DATA_DIR = config.DATA_DIR
FIGURES_DIR = config.FIGURES_DIR
OUTPUT_DATA_FILE = config.OUTPUT_DATA_FILE


def plot_output(cb_on=True):
    if not os.path.exists(FIGURES_DIR):
        os.makedirs(FIGURES_DIR)

    data = data_processing.load_sim_data()
    res = data['results']

    t = res['t']
    x1 = res['x1']
    x2 = res['x2']
    if cb_on:
        wpre = res['syn_wpre'][0]
        u = res['u'][0]
        Ca = res['Ca'][0]
        plph.plot_wpre(t, x1, wpre, u, Ca)

    # Retrieve parameters from saved metadata
    saved_params = data['params']
    num_cells = saved_params.get('NUM_CELLS', params.NUM_CELLS)
    # pph.plot_power_spec(x1, x2)
    # Generate spike matrices using loaded spike data
    spike_matrix_1 = data_processing.create_spike_matrix_histo(res['spikes_n1'], num_cells)
    spike_matrix_2 = data_processing.create_spike_matrix_histo(res['spikes_n2'], num_cells)

    pph.standard_plot(t, x1, x2, spike_matrix_1, spike_matrix_2, num_cells, params.SIM_DURATION/second+params.TRANSIENT,
                    g_inter_vals=params.G_INTER_VALS, g_intra_vals=params.G_INTRA_VALS, coupling_vals=params.COUPLING_VALS, x_naught_vals=params.X_NAUGHT_VALS)


def plot_output_full(cb_on=True):
    if not os.path.exists(FIGURES_DIR):
        os.makedirs(FIGURES_DIR)

    data = data_processing.load_sim_data()
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
        plph.plot_wpre(t, x1, wpre, u, Ca)
    # Retrieve parameters from saved metadata
    saved_params = data['params']
    num_cells = saved_params.get('NUM_CELLS', params.NUM_CELLS)

    # Generate spike matrices using loaded spike data
    spike_matrix_1 = data_processing.create_spike_matrix_histo(res['spikes_n1'], num_cells)
    spike_matrix_2 = data_processing.create_spike_matrix_histo(res['spikes_n2'], num_cells)

    pph.plot_hr_single(t, x1, y1, z1, I_syn_inter)
    pph.plot_ml_single(t, x2, n)
    pph.standard_plot(t, x1, x2, spike_matrix_1, spike_matrix_2, num_cells, params.SIM_DURATION/second+params.TRANSIENT,
                    g_inter_vals=params.G_INTER_VALS, g_intra_vals=params.G_INTRA_VALS, coupling_vals=params.COUPLING_VALS, x_naught_vals=params.X_NAUGHT_VALS)

    print("I_syn_inter max (raw amps):", np.max((I_syn_inter[0]/amp)))
    print("I_syn_intra max (raw amps):", np.max((I_syn_intra[0]/amp)))

    print("I_syn_inter min (raw amps):", np.min((I_syn_inter[0]/amp)))
    print("I_syn_intra min (raw amps):", np.min((I_syn_intra[0]/amp)))


def analyze_populations():
    data = data_processing.load_sim_data()
    res = data['results']
    x1 = res['x1']
    x2 = res['x2']
    pop1_spike_times, pop1_neuron_idx = res['spikes_n1']['t'], res['spikes_n1']['i']
    pop2_spike_times, pop2_neuron_idx = res['spikes_n2']['t'], res['spikes_n2']['i']

    data_processing.dump_spikes_to_file(np.asarray(pop1_neuron_idx), np.asarray(pop1_spike_times))

    print("============HINDMARSH ROSE STATS============")
    chi, autocorr, lag = syn.autocorrelate(x1)
    print(f'synchrony measure: {chi}\nautocorrelation: {autocorr}')
    z, r, psi = syn.KOP(pop1_neuron_idx, pop1_spike_times, params.SIM_DURATION/second)
    print(f'r: {np.mean(r)}')
    aph.plot_autocorr(autocorr, lag)
    phase_matrix = syn.compute_phase(pop1_neuron_idx, pop1_spike_times, params.SIM_DURATION/second)
    aph.plot_kop(phase_matrix)

    print("\n============MORRIS LECAR STATS============")
    chi, autocorr, lag = syn.autocorrelate(x2)
    print(f'synchrony measure: {chi}\nautocorrelation: {autocorr}')
    z, r, psi = syn.KOP(pop2_neuron_idx, pop2_spike_times, params.SIM_DURATION/second)
    print(f'r: {np.mean(r)}')


def test():
    data = data_processing.load_sim_data()
    res = data['results']
    x1 = res['x1']
    t = res['t']
    pop1_spike_times, pop1_neuron_idx = res['spikes_n1']['t'], res['spikes_n1']['i']

    pph.plot_hr_multiple(t, x1, zoom=True)
    aph.plot_auto_lfp(x1)
    phase_matrix = syn.compute_phase(pop1_neuron_idx, pop1_spike_times, params.SIM_DURATION/second)
    aph.plot_kop(phase_matrix)


def main():
    parser = argparse.ArgumentParser(description="Run and/or plot the simulation.")
    parser.add_argument('-m', '--mode', type=str, default='rp',
                        help="Run mode: 'r' run, 'p' plot, 'a' analyze, 't' test.")
    parser.add_argument('--cb', action='store_true', default=True,
                        help="Enable CB synapses (default: enabled)")
    parser.add_argument('--no-cb', dest='cb', action='store_false',
                        help="Disable CB synapses")
    args = parser.parse_args()
    run_mode = args.mode
    cb_on = args.cb

    if ('r' in run_mode):
        print("Running simulation...")
        run_sim(cb_on)
        print("Simulation complete.")
    if ('p' in run_mode):
        print("Generating plots...")
        if 'f' in run_mode:
            plot_output_full(cb_on)
        else:
            plot_output(cb_on)
        print(f"Plots saved to 'figures' directory.")
    if ('a' in run_mode):
        analyze_populations()
    if ('t' in run_mode):
        test()

if __name__ == "__main__":
    main()
