from brian2 import *
from brian2tools import *
from matplotlib import pyplot as plt
import numpy as np
import scipy
from scipy import signal
import argparse
import os
import plotting as ph # (Plotting help idk what to call it)
import config

DATA_DIR = config.DATA_DIR
FIGURES_DIR = config.FIGURES_DIR
OUTPUT_DATA_FILE = config.OUTPUT_DATA_FILE

def run_sim():
    sim_duration = 60 * second
    # change defaultclock to fix raster smearing
    tau = 1 * msecond
    defaultclock.dt = tau / 20
    print("defaultclock.dt is: ", defaultclock.dt)
    num_cells = 40

    # ISOLATE = 0: populations are decoupled (independent)
    # ISOLATE = 1: populations are coupled (normal mode)
    ISOLATE = 1

    # Discrepancy in paper regarding Wmax
    # Table shows [1, 20] with no units, while text states [2.5, 20] * mV
    Wmax = 0.006

    coupling = 1

    ### Population 1 (Hindmarsh-Rose)

    # Population 1 parameters - Set for chaotic bursting behavior
    # Classic Hindmarsh-Rose parameters for independent operation
    a = 1.0
    b = 3.0
    c = 1.0
    d = 5.0
    s = 8.0  # codebase default
    I_app_1 = 3.1  # Applied current in bursting regime (3.0-3.5 for bursting)
    x_naught = -2  # Resting potential
    r = 0.0002 / msecond  # Slow adaptation (typical: 0.001-0.01)
    sigma_1 = 1/50 # No noise coupling when isolated
    
    # Population 1 equations
    # When ISOLATE=0: all coupling (intra and inter) is removed
    # When ISOLATE=1: normal coupled behavior

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

    # refractory based on voltage
    N1 = NeuronGroup(num_cells, pop1_eqs, method='euler', threshold='x > 1.5', reset='')
    # Randomize initial values
    N1.x = np.ones(num_cells) * x_naught + randn(num_cells) * Wmax
    N1.y = 'c - d*x**2'
    N1.z = '(s*(x - x_naught))'

    # Population 1 state variable averaging
    exc_averaging_eqs ='''
    x_bar_post = x_pre / num_cells : 1 (summed)
    z_bar_post = z_pre / num_cells : 1 (summed)
    '''

    gap_junctions_1 = Synapses(N1, N1, exc_averaging_eqs)
    gap_junctions_1.connect()

    ### Population 2 (Morris-Lecar)
    # Population 2 parameters - Set for oscillatory/repetitive spiking behavior
    # Classic Morris-Lecar Type I (SNIC) parameters
    Cm = 20 * ufarad  # Membrane capacitance
    I_app_2 = 40 * uamp  # Applied current for sustained spiking -- 80 for sepeparate morris lecar
    v1 = -1.2 * mvolt  # Ca activation midpoint
    v2 = 18 * mvolt  # Ca activation slope
    v3 = 12 * mvolt  # K activation midpoint
    v4 = 17.4 * mvolt  # K activation slope (larger = slower activation)
    phi = 0.067 / msecond  # K channel time scale (faster for more regular spiking)
    E_Ca = 120 * mvolt  # Ca reversal potential
    E_K = -84 * mvolt  # K reversal potential
    E_L = -60 * mvolt  # Leak reversal potential
    gL = 2 * msiemens  # Leak conductance
    gCa = 4.0 * msiemens  # Ca conductance
    gK = 8.0 * msiemens  # K conductance
    sigma_2 = 50 * uA  
    
    # Population 2 equations    
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
    
    N2 = NeuronGroup(num_cells, pop2_eqs, method='euler', threshold='x > 0.95', reset='')
    N2.v = E_L * np.ones(num_cells) + randn(num_cells) * Wmax * volt
    N2.n = 'n_inf'
    
    # Population 2 state variable averaging
    inh_averaging_eqs ='''
    x_bar_post = x_pre / num_cells : 1 (summed)
    '''

    gap_junctions_2 = Synapses(N2, N2, inh_averaging_eqs)
    gap_junctions_2.connect()

    hindmarsh_rose_syn_eqs ='''
    du/dt = (alpha * T * (1 - u) - beta * u) : 1 (clock-driven)
    T = Tmax / (1 + exp(-(x_bar_pre * volt - Vt) / Kp)) : mM

    G : siemens
    E : volt
    alpha : mmolar ** -1 * second ** -1
    beta : second ** -1
    '''

    morris_lecar_syn_eqs='''
    du/dt = (alpha * T * (1 - u) - beta * u) : 1 (clock-driven)
    T = Tmax / (1 + exp(-(x_bar_pre * volt - Vt) / Kp)) : mM

    G : siemens
    E : volt
    alpha : mmolar ** -1 * second ** -1
    beta : second ** -1
    '''

    hr_intra_syn_eqs = '''
    I_syn_intra_post = (-G * u * (x_post * volt - E)) : amp (summed)
    ''' + hindmarsh_rose_syn_eqs

    hr_inter_syn_eqs = '''
    I_syn_inter_post = (-G * u * (x_post * volt - E)) : amp (summed)
    ''' + hindmarsh_rose_syn_eqs

    ml_intra_syn_eqs = '''
    I_syn_intra_post = (-G * u * (x_post * volt - E)) : amp (summed)
    ''' +  morris_lecar_syn_eqs

    ml_inter_syn_eqs = '''
    I_syn_inter_post = (-G * u * (x_post * volt - E)) : amp (summed)
    ''' + morris_lecar_syn_eqs

    # Synapse parameters
    Vt = 2 * mV
    Kp = 5 * mV
    Tmax = 1 * mmolar
    alpha_exc = 1.1 / (mmolar * msecond)
    alpha_inh = 5 / (mmolar * msecond)
    beta_exc = 0.19 / msecond
    beta_inh = 0.18 / msecond
    Esyn_exc = 0 * mV
    Esyn_inh = -80 * mV
    G_intra = 0.1 * uS
    G_inter = 0.2 * uS

    # Population 1 synapses to self
    S1_to_1 = Synapses(N1, N1, hr_intra_syn_eqs, method='euler')
    S1_to_1.connect()
    S1_to_1.E = Esyn_exc
    S1_to_1.alpha = alpha_exc
    S1_to_1.beta = beta_exc
    S1_to_1.G = G_intra

    # Population 1 synapses to pop 2
    S1_to_2 = Synapses(N1, N2, hr_inter_syn_eqs, method='euler')
    S1_to_2.connect()
    S1_to_2.run_regularly('z_bar_post = z_bar_pre', dt=defaultclock.dt)
    S1_to_2.E = Esyn_exc
    S1_to_2.alpha = alpha_exc
    S1_to_2.beta = beta_exc
    S1_to_2.G = G_inter

    # Population 2 synapses to self
    S2_to_2 = Synapses(N2, N2, ml_intra_syn_eqs, method='euler')
    S2_to_2.connect()
    S2_to_2.E = Esyn_inh
    S2_to_2.alpha = alpha_inh
    S2_to_2.beta = beta_inh
    S2_to_2.G = G_intra

    # Population 2 synapses to pop 1
    S2_to_1 = Synapses(N2, N1, ml_inter_syn_eqs, method='euler')
    S2_to_1.connect()
    S2_to_1.run_regularly('x2_bar_post = x_bar_pre', dt=defaultclock.dt)
    S2_to_1.E = Esyn_inh
    S2_to_1.alpha = alpha_inh
    S2_to_1.beta = beta_inh
    S2_to_1.G = G_inter

    # run(1 * second)

    # Neuron group state monitors
    # Can set dt param to record at a different time step
    M_N1 = StateMonitor(N1, ['x', 'y', 'z', 'I_syn_inter'], record=True)
    M_N2 = StateMonitor(N2, ['x', 'n', 'I_syn_inter'], record=True)

    SM_N1 = SpikeMonitor(N1)
    SM_N2 = SpikeMonitor(N2)
    
    run(sim_duration)
    
    
    
    t = np.asarray(M_N1.t)
    x1 = np.asarray(M_N1.x)
    y1 = np.asarray(M_N1.y)
    z1 = np.asarray(M_N1.z)
    I_syn_inter_1 = np.asarray(M_N1.I_syn_inter)
    x2 = np.asarray(M_N2.x)
    n2 = np.asarray(M_N2.n)
    I_syn_inter_2 = np.asarray(M_N2.I_syn_inter)

    # Save output data
    ph.save_data(OUTPUT_DATA_FILE, t=t, x1=x1, y1=y1, z1=z1, I_syn_inter_1=I_syn_inter_1, x2=x2, n2=n2, I_syn_inter_2=I_syn_inter_2)
    ph.save_data("Spike_Monitor_N1.npz", t=SM_N1.t, i=SM_N1.i)
    ph.save_data("Spike_Monitor_N2.npz", t=SM_N2.t, i=SM_N2.i)


    

def plot_output():
    if not os.path.exists(FIGURES_DIR):
        os.makedirs(FIGURES_DIR)
    
    arrs = np.load(os.path.join(DATA_DIR, OUTPUT_DATA_FILE))

    t = arrs['t']
    x1 = arrs['x1']
    y1 = arrs['y1']
    z1 = arrs['z1']
    I_syn_inter_1 = arrs['I_syn_inter_1']
    x2 = arrs['x2']
    n2 = arrs['n2']
    I_syn_inter_2 = arrs['I_syn_inter_2']

    ph.plot_both(t, x1, x2)
    ph.plot_both_avg(t, x1, y1, z1, x2, n2)
    ph.plot_hr_single(t, x1, y1, z1, I_syn_inter_1)
    ph.plot_hr_mean(t, x1, y1, z1)
    ph.plot_ml_single(t, x2, n2)

    spike_matrix_1 = create_spike_matrix_histo("Spike_Monitor_N1.npz")
    spike_matrix_2 = create_spike_matrix_histo("Spike_Monitor_N2.npz")
    ph.pop1("pop1", t, x1, spike_matrix_1)
    ph.pop1("pop2", t, x2, spike_matrix_2)

    # ph.plot_raster("N1", "Spike_Monitor_N1.npz", t, x1)
    # ph.plot_raster("N2", "Spike_Monitor_N2.npz", t, x2)

def eda():
    arrs = np.load(os.path.join(DATA_DIR, OUTPUT_DATA_FILE))

    t = arrs['t']
    x1 = arrs['x1']
    y1 = arrs['y1']
    z1 = arrs['z1']
    I_syn_inter_1 = arrs['I_syn_inter_1']
    x2 = arrs['x2']
    n2 = arrs['n2']
    I_syn_inter_2 = arrs['I_syn_inter_2']

    print(f"Length of t is: {t.shape}")
    print(f"Length of x1 is: {x1.shape}")
    print(f"Length of ISYNINTER is: {I_syn_inter_1.shape}")
    print(f"Length of x2 is: {x2.shape}")

    bin_freq = 100
    num_ticks = bin_freq * 20
    meow = len(t)//num_ticks
    new_x1 = np.array([x1[:, i*num_ticks:i*num_ticks+num_ticks].mean(axis=1) for i in range(meow)]).T
    new_x2 = np.array([x2[:, i*num_ticks:i*num_ticks+num_ticks].mean(axis=1) for i in range(meow)]).T
    new_t = np.array([t[i*num_ticks] for i in range(meow)]).T
    # new_x1 = x1.reshape(x1.shape[0], -1, num_ticks).mean(axis=2)
    print(new_t.shape)

# this does the same as the histogram version but worse
def create_spike_matrix(data_name):
    arrs = np.load(os.path.join(DATA_DIR, data_name))
    spike_times = arrs['t']  # All spike times
    neuron_indices = arrs['i']  # Corresponding neuron indices
    # Define time bins
    sim_duration = 60 
    dt = 0.1  # 100ms
    time_bins = np.arange(0, sim_duration, dt)
    num_bins = len(time_bins)
    num_neurons = 40
    
    # Create empty spike matrix
    spike_matrix = np.zeros((num_neurons, num_bins))
    
    # Bin the spikes
    bin_indices = np.digitize(spike_times, time_bins)-1
    
    # Fill the matrix
    for idx in range(len(bin_indices)):
        bin_idx = bin_indices[idx]
        if bin_idx <= 15:
            # filter out initial warmup 'spikes'
            continue
        neuron = neuron_indices[idx]
        spike_matrix[neuron][bin_idx] += 1

    return spike_matrix

def create_spike_matrix_histo(data_name):
    """Create spike matrix for imshow to plot raster plots"""
    arrs = np.load(os.path.join(DATA_DIR, data_name))
    spike_times = arrs['t']  # All spike times
    neuron_indices = arrs['i']  # Corresponding neuron indices

    # Define time bins
    sim_duration = 60
    dt = 0.01  # 10ms
    warmup_time = 2  # Skip warup period - spikes dont count

    # Filter warmup spikes
    valid = spike_times > warmup_time
    spike_times = spike_times[valid]
    neuron_indices = neuron_indices[valid]

    num_neurons = 40

    # Create bin edges (need +1 for right edge)
    time_bins = np.arange(0, sim_duration + dt, dt)
    neuron_bins = np.arange(0, num_neurons + 1)

    # x = rows, y = columns -- this is flipped from math coordinates, so swap everything accordingly 
    spike_matrix, neuron_edges, time_edges = np.histogram2d(
        neuron_indices, # x
        spike_times,    # y
        bins=[neuron_bins, time_bins]
    )

    return spike_matrix


def main():
    ### Run mode string
    # r - run simulation
    # p - plot results
    parser = argparse.ArgumentParser(description="Run and/or plot the simulation.")
    parser.add_argument('-m', '--mode', type=str, default='rp', 
                        help="Run mode: 'r' to run, 'p' to plot, 'rp' to run and plot. Default is 'rp'.")
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
    if ('t' in run_mode):
        # testing features
        create_spike_matrix_histo()

if __name__ == "__main__":
    main()
