from brian2 import *
from brian2tools import *
import numpy as np
import argparse
import os
import config
import params 
import plotting as ph 
import synch as syn
import data_processing

DATA_DIR = config.DATA_DIR
FIGURES_DIR = config.FIGURES_DIR
OUTPUT_DATA_FILE = config.OUTPUT_DATA_FILE

def run_sim():
    # Setup Simulation
    defaultclock.dt = params.TAU_CLOCK / params.DT_SCALING
    print("defaultclock.dt is: ", defaultclock.dt)
    
    # --- Population 1: Hindmarsh-Rose ---
    pop1_eqs = '''
    dx/dt = (y - a * x ** 3 + b * x ** 2 - z + I_app_1
        + ISOLATE * (coupling * (x_bar - x)
        + Wmax * xi * sqrt(second)
        + sigma_1 * (I_syn_intra + I_syn_inter) / amp)) / tau : 1
    dy/dt = (c - d * x ** 2 - y) / tau : 1
    dz/dt = r * (s * (x + ISOLATE * x2_bar - x_naught) - z_bar) : 1

    x_bar : 1
    z_bar : 1
    x2_bar : 1
    I_syn_intra : amp
    I_syn_inter : amp
    '''

    pop1_namespace = {
        'a': params.HR_A, 'b': params.HR_B, 'c': params.HR_C, 
        'd': params.HR_D, 's': params.HR_S, 'I_app_1': params.HR_I_APP,
        'x_naught': params.HR_X_NAUGHT, 'r': params.HR_R, 
        'sigma_1': params.HR_SIGMA, 'tau': params.TAU_CLOCK,
        'ISOLATE': params.ISOLATE, 'coupling': params.COUPLING_STRENGTH,
        'Wmax': params.W_MAX
    }

    N1 = NeuronGroup(params.NUM_CELLS, pop1_eqs, method='euler', 
                     threshold=params.HR_THRESHOLD, reset='',
                     namespace=pop1_namespace, refractory=params.HR_REFRACTORY_CONDITION)
    
    N1.x = np.ones(params.NUM_CELLS) * params.HR_X_NAUGHT + randn(params.NUM_CELLS) * params.W_MAX
    N1.y = 'c - d*x**2'
    N1.z = '(s*(x - x_naught))'

    # Population 1 Averaging
    exc_averaging_eqs ='''
    x_bar_post = x_pre / num_cells : 1 (summed)
    z_bar_post = z_pre / num_cells : 1 (summed)
    '''
    gap_junctions_1 = Synapses(N1, N1, exc_averaging_eqs, namespace={'num_cells': params.NUM_CELLS})
    gap_junctions_1.connect()

    
    # --- Population 2: Morris-Lecar ---    
    pop2_eqs = '''
    dv/dt = (I_app_2 - gL*(v-E_L) - gK*n*(v-E_K) - gCa*m_inf*(v-E_Ca)
        + ISOLATE * (sigma_2 * (Wmax * xi * sqrt(second)
        + coupling * (x_bar - x) - 0.3 * (z_bar - 3))
        + (I_syn_intra + I_syn_inter))) / Cm : volt
    dn/dt = phi * (n_inf - n) / tau_n : 1

    m_inf = 0.5 * (1 + tanh((v - v1) / v2)) : 1
    n_inf = 0.5 * (1 + tanh((v - v3) / v4)) : 1
    tau_n = 1 / cosh((v - v3) / (2 * v4)) : 1
    
    x = v/(20 * mV) : 1
    
    x_bar : 1
    z_bar : 1
    I_syn_inter : amp
    I_syn_intra : amp
    '''
    
    pop2_namespace = {
        'Cm': params.ML_CM, 'I_app_2': params.ML_I_APP, 'gL': params.ML_GL,
        'E_L': params.ML_E_L, 'gK': params.ML_GK, 'E_K': params.ML_E_K,
        'gCa': params.ML_GCA, 'E_Ca': params.ML_E_CA, 
        'v1': params.ML_V1, 'v2': params.ML_V2, 'v3': params.ML_V3, 'v4': params.ML_V4,
        'phi': params.ML_PHI, 'sigma_2': params.ML_SIGMA,
        'Wmax': params.W_MAX, 'coupling': params.COUPLING_STRENGTH,
        'ISOLATE': params.ISOLATE
    }

    N2 = NeuronGroup(params.NUM_CELLS, pop2_eqs, method='euler', 
                     threshold=params.ML_THRESHOLD, reset='',
                     namespace=pop2_namespace, refractory=params.ML_REFRACTORY_CONDITION)
    
    N2.v = params.ML_E_L * np.ones(params.NUM_CELLS) + \
           randn(params.NUM_CELLS) * params.W_MAX * volt
    N2.n = 'n_inf'
    
    # Population 2 Averaging
    inh_averaging_eqs ='''
    x_bar_post = x_pre / num_cells : 1 (summed)
    '''
    gap_junctions_2 = Synapses(N2, N2, inh_averaging_eqs, namespace={'num_cells': params.NUM_CELLS})
    gap_junctions_2.connect()


    # --- Synapses ---
    syn_namespace = {
        'Tmax': params.SYN_TMAX,
        'Vt': params.SYN_VT,
        'Kp': params.SYN_KP
    }

    syn_input_scale = 1/pop1_namespace['sigma_1']
    syn_eqs ='''
    du/dt = (alpha * T * (1 - u) - beta * u) : 1 (clock-driven)
    T = Tmax / (1 + exp(-(x_bar_pre * (syn_input_scale) * mvolt - Vt) / Kp)) : mM

    G : siemens
    E : volt
    alpha : mmolar ** -1 * second ** -1
    beta : second ** -1
    '''

    intra_syn_eqs = '''
    I_syn_intra_post = (-G * u * (x_post * (syn_input_scale) * mvolt - E)) : amp (summed)
    ''' + syn_eqs

    inter_syn_eqs = '''
    I_syn_inter_post = (-G * u * (x_post * (syn_input_scale) * mvolt - E)) : amp (summed)
    ''' + syn_eqs

    S1_to_1 = Synapses(N1, N1, intra_syn_eqs, method='euler', namespace=syn_namespace)
    S1_to_1.connect()
    S1_to_1.E = params.SYN_E_EXC
    S1_to_1.alpha = params.SYN_ALPHA_EXC
    S1_to_1.beta = params.SYN_BETA_EXC
    S1_to_1.G = params.G_INTRA

    S1_to_2 = Synapses(N1, N2, inter_syn_eqs, method='euler', namespace=syn_namespace)
    S1_to_2.connect()
    S1_to_2.run_regularly('z_bar_post = z_bar_pre', dt=defaultclock.dt)
    S1_to_2.E = params.SYN_E_EXC
    S1_to_2.alpha = params.SYN_ALPHA_EXC
    S1_to_2.beta = params.SYN_BETA_EXC
    S1_to_2.G = params.G_INTER

    S2_to_2 = Synapses(N2, N2, intra_syn_eqs, method='euler', namespace=syn_namespace)
    S2_to_2.connect()
    S2_to_2.E = params.SYN_E_INH
    S2_to_2.alpha = params.SYN_ALPHA_INH
    S2_to_2.beta = params.SYN_BETA_INH
    S2_to_2.G = params.G_INTRA

    S2_to_1 = Synapses(N2, N1, inter_syn_eqs, method='euler', namespace=syn_namespace)
    S2_to_1.connect()
    S2_to_1.run_regularly('x2_bar_post = x_bar_pre', dt=defaultclock.dt)
    S2_to_1.E = params.SYN_E_INH
    S2_to_1.alpha = params.SYN_ALPHA_INH
    S2_to_1.beta = params.SYN_BETA_INH
    S2_to_1.G = params.G_INTER

    M_N1 = StateMonitor(N1, ['x', 'y', 'z', 'I_syn_inter'], record=True)
    M_N2 = StateMonitor(N2, ['x', 'n', 'I_syn_inter'], record=True)

    SM_N1 = SpikeMonitor(N1)
    SM_N2 = SpikeMonitor(N2)
    
    
    # Run
    run(params.SIM_DURATION)
    data_processing.save_data(M_N1, M_N2, SM_N1, SM_N2)

def plot_output():
    if not os.path.exists(FIGURES_DIR):
        os.makedirs(FIGURES_DIR)
    
    data = data_processing.load_sim_data()
    res = data['results']
    t = data_processing.cutoff_transient(res['t'], params.TRANSIENT, params.TAU_CLOCK/params.DT_SCALING/msecond*1e-3)
    x1 = data_processing.cutoff_transient(res['x1'],  params.TRANSIENT, params.TAU_CLOCK/params.DT_SCALING/msecond*1e-3)
    x2 = data_processing.cutoff_transient(res['x2'],  params.TRANSIENT, params.TAU_CLOCK/params.DT_SCALING/msecond*1e-3)

    # Retrieve parameters from saved metadata
    saved_params = data['params']
    num_cells = saved_params.get('NUM_CELLS', params.NUM_CELLS)

    # Generate spike matrices using loaded spike data
    spike_matrix_1 = data_processing.create_spike_matrix_histo(res['spikes_n1'], num_cells,  params.TRANSIENT)
    spike_matrix_2 = data_processing.create_spike_matrix_histo(res['spikes_n2'], num_cells,  params.TRANSIENT)

    ph.standard_plot(t, x1, x2, spike_matrix_1, spike_matrix_2, num_cells, params.SIM_DURATION/second)


def plot_output_full():
    if not os.path.exists(FIGURES_DIR):
        os.makedirs(FIGURES_DIR)
    
    data = data_processing.load_sim_data()
    res = data['results']
    t = res['t']
    x1 = res['x1']
    y1 = res['y1']
    z1 = res['z1']
    I_syn_inter_1 = res['I_syn_inter_1']
    x2 = res['x2']
    n = res['n2']
    
    # Retrieve parameters from saved metadata
    saved_params = data['params']
    num_cells = saved_params.get('NUM_CELLS', params.NUM_CELLS)
    
    # Generate spike matrices using loaded spike data
    spike_matrix_1 = data_processing.create_spike_matrix_histo(res['spikes_n1'], num_cells, params.TRANSIENT)
    spike_matrix_2 = data_processing.create_spike_matrix_histo(res['spikes_n2'], num_cells, params.TRANSIENT)

    ph.plot_hr_single(t, x1, y1, z1, I_syn_inter_1)
    ph.plot_ml_single(t, x2, n)
    ph.standard_plot(t, x1, x2, spike_matrix_1, spike_matrix_2, num_cells, params.SIM_DURATION/second)


def analyze_populations():
    data = data_processing.load_sim_data()
    res = data['results']
    x1 = res['x1']
    x2 = res['x2']
    pop1_spike_data = res['spikes_n1']
    pop1_spike_times, pop1_neuron_idx = pop1_spike_data['t'], pop1_spike_data['i']
    pop2_spike_data = res['spikes_n2']
    pop2_spike_times, pop2_neuron_idx = pop2_spike_data['t'], pop2_spike_data['i']

    data_processing.dump_spikes_to_file(np.asarray(pop1_neuron_idx), np.asarray(pop1_spike_times))


    print("============HINDMARSH ROSE STATS============")
    chi, autocorr = syn.autocorelate(x1)
    print(f'synchrony measure: {chi}\nautocorrelation: {autocorr}')
    z, r, psi = syn.KOP(pop1_neuron_idx, pop1_spike_times, params.SIM_DURATION/second)
    print(f'z: {z}')
    print(f'r: {r}')
    print(f'psi: {psi}')

    print("\n============MORRIS LECAR STATS============")
    chi, autocorr = syn.autocorelate(x2)
    print(f'synchrony measure: {chi}\nautocorrelation: {autocorr}')
    z, r, psi = syn.KOP(pop2_neuron_idx, pop2_spike_times, params.SIM_DURATION/second)
    print(f'z: {z}')
    print(f'r: {r}')
    print(f'psi: {psi}')



def main():
    parser = argparse.ArgumentParser(description="Run and/or plot the simulation.")
    parser.add_argument('-m', '--mode', type=str, default='rp', 
                        help="Run mode: 'r' to run, 'p' to plot, 'rp' to run and plot.")
    args = parser.parse_args()
    run_mode = args.mode

    if ('r' in run_mode):
        print("Running simulation...")
        run_sim()
        print("Simulation complete.")
    if ('p' in run_mode):
        print("Generating plots...")
        if 'f' in run_mode:
            plot_output_full()
        else:
            plot_output()
        print(f"Plots saved to 'figures' directory.")
    if ('s' in run_mode):
        analyze_populations()

if __name__ == "__main__":
    main()