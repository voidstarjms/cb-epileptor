from brian2 import *
from brian2tools import *
from matplotlib import pyplot as plt
import numpy as np
import scipy
from scipy import signal
import argparse
import os

def run_sim():
    sim_duration = 15 * second
    tau = 1 * msecond
    defaultclock.dt = tau / 20

    num_cells = 40

    # Discrepancy in paper regarding Wmax
    # Table shows [1, 20] with no units, while text states [2.5, 20] * mV
    Wmax = 0.0025

    coupling = 1

    ### Population 1 (Hindmarsh-Rose)

    # Population 1 parameters
    # NOTE: these parameters are from Brian docs, not Naze paper
    a = 1
    b = 3
    c = 1
    d = 5
    s = 4
    I_app_1 = 2
    x_naught = -1.6
    r = 0.001
    sigma_1 = 1/50
    
    # Population 1 equations
    pop1_eqs = '''
    dx/dt = (y - a * x ** 3 + b * x ** 2 - z + I_app_1
        + sigma_1 * (coupling * (x_bar - x)
        + 2 * Wmax * xi * sqrt(second) - Wmax
        + (I_syn_intra + I_syn_inter) / amp)) / tau : 1
    dy/dt = (c - d * x ** 2 - y) / tau : 1
    dz/dt = r * (s * (x + x2_bar - x_naught) - z_bar) / tau : 1

    x_bar : 1
    z_bar : 1
    x2_bar : 1
    I_syn_intra : amp
    I_syn_inter : amp
    '''

    N1 = NeuronGroup(num_cells, pop1_eqs, method='euler', threshold='x > 1.9', reset='')
    # Randomize initial values
    N1.x = np.ones(num_cells) * x_naught + randn(num_cells) * Wmax
    N1.y = 'c - d*x**2'
    N1.z = 'r*(s*(x - x_naught))'

    # Population 1 gap junctions
    gap_junction_exc_eqs ='''
    x_bar_post = x_pre / num_cells : 1 (summed)
    z_bar_post = z_pre / num_cells : 1 (summed)
    '''

    gap_junctions_1 = Synapses(N1, N1, gap_junction_exc_eqs)
    gap_junctions_1.connect()

    ### Population 2 (Morris-Lecar)
    # NOTE: these parameters are from Brian docs, not Naze paper
    # Population 2 parameters
    Cm = 20 * ufarad
    I_app_2 = 400 * uamp
    v1 = 10 * mvolt
    v2 = 15 * mvolt
    v3 = -1 * mvolt
    v4 = 14.5 * mvolt
    phi = 1.0 / (15*ms)
    E_Ca = 100 * mvolt
    E_K = -70 * mvolt
    E_L = -50 * mvolt
    gL = 2 * msiemens
    gCa = 4 * msiemens
    gK = 8 * msiemens
    sigma_2 = 50 * uA
    
    # Population 2 equations    
    pop2_eqs = '''
    dv/dt = (I_app_2 - gL*(v-E_L) - gK*n*(v-E_K) - gCa*m_inf*(v-E_Ca)
        + sigma_2 * (2 * Wmax * xi * sqrt(second) - Wmax +
        coupling * (x_bar - x) - 0.3 * (z_bar - 3)) + I_syn_intra + I_syn_inter) / Cm : volt
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
    
    # Population 2 gap junctions
    gap_junction_inh_eqs ='''
    x_bar_post = x_pre / num_cells : 1 (summed)
    '''

    gap_junctions_2 = Synapses(N2, N2, gap_junction_inh_eqs)
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

    run(60 * second)

    # Neuron group state monitors
    M_N1 = StateMonitor(N1, ['x', 'y', 'z'], record=True)
    M_N2 = StateMonitor(N2, ['x', 'v'], record=True)

    SM_N1 = SpikeMonitor(N1)
    SM_N2 = SpikeMonitor(N2)
    
    run(sim_duration)
    
    t = np.asarray(M_N1.t)
    x1 = np.asarray(M_N1.x)
    y1 = np.asarray(M_N1.y)
    z1 = np.asarray(M_N1.z)
    x2 = np.asarray(M_N2.x)
    v2 = np.asarray(M_N2.v)
    #n1_spikes = np.asarray(SM_N1.spike_trains())

    brian_plot(SM_N1)
    plt.show()
    brian_plot(SM_N2)
    plt.show()

    np.savez("output_data.npz", t=t, x1=x1, y1=y1, z1=z1, x2=x2, v2=v2)

def plot_output():
    arrs = np.load("output_data.npz")

    t = arrs['t']
    x1 = arrs['x1']
    y1 = arrs['y1']
    z1 = arrs['z1']
    x2 = arrs['x2']
    v2 = arrs['v2']

    # fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    # ax1.plot(t, x1[0])
    # ax1.set_xlabel("Time (s)")
    # ax1.set_ylabel("x1")
    # ax2.plot(t, x2[0])
    # ax2.set_xlabel("Time (s)")
    # ax2.set_ylabel("x2")

    pop1_mean = np.mean(x1, axis=0)
    pop2_mean = np.mean(x2, axis=0)
    mean_potential = 0.8 * pop1_mean + 0.2 * pop2_mean
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 8))
    ax1.plot(t, pop2_mean)
    
    fs = 1 / defaultclock.dt / Hz
    f, Pxx = signal.welch(mean_potential, fs=fs)
    ax2.semilogy(f, Pxx)
    
    # Output directory
    output_dir = "figures"
    os.makedirs(output_dir, exist_ok=True)
    
    plt.savefig(os.path.join(output_dir, "interictal_power_spectrum.png"), format="png")
    plt.show()

    fig = plt.figure()
    f, ts, Sxx = scipy.signal.spectrogram(mean_potential, fs)
    fig = plt.pcolormesh(ts, f, Sxx, shading='gouraud')

    plt.savefig(os.path.join(output_dir, "interictal_spectrogram.png"), format="png")
    plt.show()

def main():
    ### Run mode string
    # r - run simulation
    # s - save simulation results
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

if __name__ == "__main__":
    main()
