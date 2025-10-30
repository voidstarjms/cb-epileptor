from brian2 import *
from brian2tools import *
from matplotlib import pyplot as plt
import numpy as np
import scipy
from scipy import signal
import argparse
import os


# File paths and directories
DATA_DIR = 'data/'
FIGURES_DIR = 'figures/'
OUTPUT_DATA_FILE = 'output_data.npz'

def run_sim():
    sim_duration = 15 * second
    tau = 1 * msecond
    defaultclock.dt = tau / 20
    print("defaultclock.dt is: ", defaultclock.dt)
    num_cells = 40

    # Discrepancy in paper regarding Wmax
    # Table shows [1, 20] with no units, while text states [2.5, 20] * mV
    Wmax = 0.004

    coupling = 1

    ### Population 1 (Hindmarsh-Rose)

    # Population 1 parameters
    # NOTE: these parameters are from Brian docs, not Naze paper
    a = 1
    b = 3
    c = 1
    d = 5
    s = 8
    I_app_1 = 3.1
    x_naught = -3
    r = 0.0003 / msecond
    sigma_1 = 1/50
    
    # Population 1 equations
    pop1_eqs = '''
    dx/dt = (y - a * x ** 3 + b * x ** 2 - z + I_app_1
        + coupling * (x_bar - x) + Wmax * xi * sqrt(second)
        + sigma_1 * (I_syn_intra + I_syn_inter) / amp) / tau : 1
    dy/dt = (c - d * x ** 2 - y) / tau : 1
    dz/dt = r * (s * (x + x2_bar - x_naught) - z_bar) : 1

    x_bar : 1
    z_bar : 1
    x2_bar : 1
    I_syn_intra : amp
    I_syn_inter : amp
    '''

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
    # NOTE: these parameters are from Brian docs, not Naze paper
    # Population 2 parameters
    Cm = 20 * ufarad
    I_app_2 = 40 * uamp
    v1 = 1.2 * mvolt
    v2 = 18 * mvolt
    v3 = 12 * mvolt
    v4 = 17.4 * mvolt
    phi = 1.0 / (15*ms)
    E_Ca = 120 * mvolt
    E_K = -84 * mvolt
    E_L = -60 * mvolt
    gL = 2 * msiemens
    gCa = 4 * msiemens
    gK = 8 * msiemens
    sigma_2 = 50 * uA
    
    # Population 2 equations    
    pop2_eqs = '''
    dv/dt = (I_app_2 - gL*(v-E_L) - gK*n*(v-E_K) - gCa*m_inf*(v-E_Ca)
        + sigma_2 * (Wmax * xi * sqrt(second)
        + coupling * (x_bar - x) - 0.3 * (z_bar - 3))
        + 20 * (I_syn_intra + I_syn_inter)) / Cm : volt
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

    syn_eqs ='''
    du/dt = (alpha * T * (1 - u) - beta * u) : 1 (clock-driven)
    T = Tmax / (1 + exp(-(x_bar_pre * mvolt - Vt) / Kp)) : mM

    G : siemens
    E : volt
    alpha : mmolar ** -1 * second ** -1
    beta : second ** -1
    '''

    intra_syn_eqs = '''
    I_syn_intra_post = (-G * u * (x_post * mvolt - E)) : amp (summed)
    ''' + syn_eqs

    inter_syn_eqs = '''
    I_syn_inter_post = (-G * u * (x_post * mvolt - E)) : amp (summed)
    ''' + syn_eqs

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
    G_intra = 0.2 * uS
    G_inter = 0.1 * uS

    # Population 1 synapses to self
    S1_to_1 = Synapses(N1, N1, intra_syn_eqs, method='euler')
    S1_to_1.connect()
    S1_to_1.E = Esyn_exc
    S1_to_1.alpha = alpha_exc
    S1_to_1.beta = beta_exc
    S1_to_1.G = G_intra

    # Population 1 synapses to pop 2
    S1_to_2 = Synapses(N1, N2, inter_syn_eqs, method='euler')
    S1_to_2.connect()
    S1_to_2.run_regularly('x2_bar_post = x_bar_pre', dt=defaultclock.dt)
    S1_to_2.run_regularly('z_bar_post = z_bar_pre', dt=defaultclock.dt)
    S1_to_2.E = Esyn_exc
    S1_to_2.alpha = alpha_exc
    S1_to_2.beta = beta_exc
    S1_to_2.G = G_inter

    # Population 2 synapses to self
    S2_to_2 = Synapses(N2, N2, intra_syn_eqs, method='euler')
    S2_to_2.connect()
    S2_to_2.E = Esyn_inh
    S2_to_2.alpha = alpha_inh
    S2_to_2.beta = beta_inh
    S2_to_2.G = G_intra

    # Population 2 synapses to pop 1
    S2_to_1 = Synapses(N2, N1, inter_syn_eqs, method='euler')
    S2_to_1.connect()
    S2_to_1.E = Esyn_inh
    S2_to_1.alpha = alpha_inh
    S2_to_1.beta = beta_inh
    S2_to_1.G = G_inter

    # run(1 * second)

    # Neuron group state monitors
    M_N1 = StateMonitor(N1, ['x', 'y', 'z'], record=True)
    M_N2 = StateMonitor(N2, ['x', 'n'], record=True)

    SM_N1 = SpikeMonitor(N1)
    SM_N2 = SpikeMonitor(N2)
    
    run(sim_duration)
    
    
    
    t = np.asarray(M_N1.t)
    x1 = np.asarray(M_N1.x)
    y1 = np.asarray(M_N1.y)
    z1 = np.asarray(M_N1.z)
    x2 = np.asarray(M_N2.x)
    n2 = np.asarray(M_N2.n)

    #n1_spikes = np.asarray(SM_N1.spike_trains())


    # Save output data
    save_data(OUTPUT_DATA_FILE, t=t, x1=x1, y1=y1, z1=z1, x2=x2, n2=n2)
    save_data("Spike_Monitor_N1.npz", t=SM_N1.t, i=SM_N1.i)
    save_data("Spike_Monitor_N2.npz", t=SM_N2.t, i=SM_N2.i)

# kwargs collects args into a dict, allows flexible arguments to be passed in
def save_data(filename, **kwargs):
    if not os.path.exists(DATA_DIR):
        os.makedirs(DATA_DIR)
    
    np.savez(os.path.join(DATA_DIR, filename), **kwargs)

def plot_raster(fig_name, data_filename):
    if not os.path.exists(FIGURES_DIR):
        os.makedirs(FIGURES_DIR)
    moni = np.load(os.path.join(DATA_DIR, data_filename))
    
    plt.figure(figsize=(12, 8))
    plt.plot(moni['t']/ms, moni['i'], '.k', markersize=2)
    plt.title(f'{fig_name} - Raster')
    plt.xlabel('Time (ms)')
    plt.ylabel('Neuron Index')
    plt.grid(True, alpha=0.3)
    plt.savefig(os.path.join(FIGURES_DIR, f"{fig_name}_raster.png"), format="png", dpi=300, bbox_inches='tight')
    plt.show()

def plot_output():
    if not os.path.exists(FIGURES_DIR):
        os.makedirs(FIGURES_DIR)
    
    arrs = np.load(os.path.join(DATA_DIR, OUTPUT_DATA_FILE))

    t = arrs['t']
    x1 = arrs['x1']
    y1 = arrs['y1']
    z1 = arrs['z1']
    x2 = arrs['x2']
    n2 = arrs['n2']

    # One neuron from both pops 
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    ax1.plot(t, x1[0])
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("x1")
    ax2.plot(t, x2[0])
    ax2.set_xlabel("Time (s)")
    ax2.set_ylabel("x2")

    # All pop 1 variables
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(30, 10), sharex=True)
    ax1.plot(t, x1[0])
    ax1.set_ylabel("Neuron 0 x")
    ax2.plot(t, y1[0])
    ax2.set_ylabel("Neuron 0 y")
    ax3.plot(t, z1[0])
    ax3.set_ylabel("Neuron 0 z")
    ax3.set_xlabel("Time (s)")

    # All pop2 variables
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(30, 10), sharex=True)
    ax1.plot(t, x2[0])
    ax1.set_ylabel("Neuron 0 x")
    ax2.plot(t, n2[0])
    ax2.set_ylabel("Neuron 0 n")
    ax2.set_xlabel("Time (s)")
    
    #plt.savefig("figures/interictal_pop2_r4e-5_10s.png", format="png")
    plt.savefig(os.path.join(FIGURES_DIR, "interictal_pop2_test.png"), format="png")
    plt.show()
    
    # pop1_mean = np.mean(x1, axis=0)
    # pop2_mean = np.mean(x2, axis=0)
    # mean_potential = 0.8 * pop1_mean + 0.2 * pop2_mean
    # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 8))
    # ax1.plot(t, pop2_mean)
    # ax1.set_xlabel("Time (s)")
    # ax1.set_ylabel("Weighted mean potential (a.u.)")
    
    # fs = 1 / defaultclock.dt / Hz
    # f, Pxx = signal.welch(mean_potential, fs=fs)
    # ax2.semilogy(f, Pxx)
    # ax2.set_xlabel("Frequency (Hz)")
    # ax2.set_ylabel("Amplitude (a.u.)")
    # plt.savefig("figures/interictal_power_spectrum_hi_r.png", format="png")
    # plt.show()

    # fig = plt.figure()
    # f, ts, Sxx = scipy.signal.spectrogram(mean_potential, fs)
    # fig = plt.pcolormesh(ts, f, Sxx, shading='gouraud')


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
        plot_raster("N1", "Spike_Monitor_N1.npz")
        plot_raster("N2", "Spike_Monitor_N2.npz")
        print(f"Plots saved to 'figures' directory.")
        

if __name__ == "__main__":
    main()
