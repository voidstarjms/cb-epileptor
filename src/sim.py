from brian2 import *
from brian2tools import *
from matplotlib import pyplot as plt
import numpy as np
import argparse
import os
import datetime 
import pickle
import plotting as ph 
import config
import params 

DATA_DIR = config.DATA_DIR
FIGURES_DIR = config.FIGURES_DIR
OUTPUT_DATA_FILE = config.OUTPUT_DATA_FILE

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
                     namespace=pop1_namespace)
    
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
                     namespace=pop2_namespace)
    
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

    syn_input_scale = 1/pop1_namespace['sigma_1'] * 20
    print(f'======={syn_input_scale}========')
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


def plot_output():
    if not os.path.exists(FIGURES_DIR):
        os.makedirs(FIGURES_DIR)
    
    # Load pickle
    filepath = os.path.join(DATA_DIR, OUTPUT_DATA_FILE)
    with open(filepath, 'rb') as f:
        data = pickle.load(f)

    # Unpack results
    res = data['results']
    t = res['t']
    x1 = res['x1']
    x2 = res['x2']
    
    # Retrieve parameters from saved metadata
    saved_params = data['params']
    num_cells = saved_params.get('NUM_CELLS', params.NUM_CELLS)
    
    # Generate spike matrices using loaded spike data
    spike_matrix_1 = create_spike_matrix_histo(res['spikes_n1'], num_cells)
    spike_matrix_2 = create_spike_matrix_histo(res['spikes_n2'], num_cells)

    ph.standard_plot(t, x1, x2, spike_matrix_1, spike_matrix_2, num_cells, params.SIM_DURATION/second)


def create_spike_matrix_histo(spike_data, num_cells):
    spike_times = spike_data['t'] 
    neuron_indices = spike_data['i'] 

    duration = params.SIM_DURATION/second
    dt = 0.1  # 100ms per bin
    warmup_time = 0

    valid = spike_times > warmup_time
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
        plot_output()
        print(f"Plots saved to 'figures' directory.")

if __name__ == "__main__":
    main()