from brian2 import *
from matplotlib import pyplot as plt
import argparse

def main():
    sim_duration = 0.2 * second
    tau = 1 * msecond
    defaultclock.dt = tau/200

    num_cells = 40

    Wmax = 0#20

    ### Population 1 (Hindmarsh-Rose)

    # Population 1 parameters
    # NOTE: these parameters are from Brian docs, not Naze paper
    a = 1
    b = 3
    c = 1
    d = 5
    s = 8
    I_app_1 = 3.1
    x_naught = -4.5
    r = 0.00002

    # Population 1 equations
    pop1_eqs = '''
    dx/dt = (y - a * x ** 3 + b * x ** 2 - z + I_app_1 + sigma_1 * (x_bar - x)
        + sigma_1 * (2 * Wmax * xi * sqrt(second) - Wmax)
        + sigma_1 * I_syn_intra / amp + sigma_1 * I_syn_inter / amp) / tau : 1
    dy/dt = (c - d * x ** 2 - y) / tau : 1
    dz/dt = r * (s * (x + x2_bar - x_naught) - z_bar) / tau : 1

    x_bar : 1
    z_bar : 1
    x2_bar : 1
    I_syn_intra : amp
    I_syn_inter : amp
    '''

    N1 = NeuronGroup(num_cells, pop1_eqs, method='euler') 
    N1.x = x_naught
    N1.y = 'c - d*x**2'
    N1.z = 'r*(s*(x - x_naught))'

    # Population 1 gap junctions
    sigma_1 = 1/50
    gap_junction_exc_eqs ='''
    x_bar_post = x / num_cells : 1 (summed)
    z_bar_post = z / num_cells : 1 (summed)
    '''

    gap_junction_inh_eqs ='''
    x_bar_post = x / num_cells : 1 (summed)
    '''

    gap_junctions_1 = Synapses(N1, N1, gap_junction_exc_eqs)
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

    # Population 2 equations    
    pop2_eqs = '''
    dv/dt = (I_app_2 - gL*(v-E_L) - gK*n*(v-E_K) - gCa*m_inf*(v-E_Ca)
        + sigma_2 * (2 * Wmax * xi * sqrt(second) - Wmax) +
        sigma_2 * (x_bar - x) + I_syn_inter + I_syn_intra
        - sigma_2 * 0.3 * (z_bar - 3)) / Cm : volt
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
    
    N2 = NeuronGroup(num_cells, pop2_eqs, method='euler')
    N2.v = E_L
    N2.n = 'n_inf'
    
    # Population 2 gap junctions
    sigma_2 = 50 * uA
    gap_junctions_2 = Synapses(N2, N2, gap_junction_inh_eqs)
    gap_junctions_2.connect()

    universal_syn_eqs ='''
    du/dt = (alpha * T * (1 - u) - beta * u) : 1 (clock-driven)
    T = Tmax / (1 + exp(-(x_pre * mvolt - Vt) / Kp)) : mM

    alpha : mmolar ** -1 * second ** -1
    beta : second ** -1
    '''

    syn_intra_eqs ='''
    I_syn_intra_post = (-G * u * (x_post * mvolt - E)) : amp (summed)

    G : siemens
    E : volt
    ''' + universal_syn_eqs

    syn_inter_inh_eqs ='''
    I_syn_inter_post = (-G * u * (x_post * mvolt - E)) : amp (summed)

    #x2_bar_post = x_bar_pre : 1

    G : siemens
    E : volt
    ''' + universal_syn_eqs

    syn_inter_exc_eqs ='''
    #z_bar_post = z_bar_pre : 1
    ''' + syn_inter_inh_eqs

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

    # Population 1 synapses to self
    S1_to_1 = Synapses(N1, N1, syn_intra_eqs)
    S1_to_1.connect()
    S1_to_1.E = Esyn_exc
    S1_to_1.alpha = alpha_exc
    S1_to_1.beta = beta_exc

    # Population 1 synapses to pop 2
    S1_to_2 = Synapses(N1, N2, syn_inter_exc_eqs)
    S1_to_2.connect()
    S1_to_2.run_regularly('x2_bar_post = x_bar_pre', dt=defaultclock.dt)
    S1_to_2.run_regularly('z_bar_post = z_bar_pre', dt=defaultclock.dt)
    S1_to_2.E = Esyn_exc
    S1_to_2.alpha = alpha_exc
    S1_to_2.beta = beta_exc

    # Population 2 synapses to self
    S2_to_2 = Synapses(N2, N2, syn_intra_eqs)
    S2_to_2.connect()
    S2_to_2.E = Esyn_inh
    S2_to_2.alpha = alpha_inh
    S2_to_2.beta = beta_inh

    # Population 2 synapses to pop 1
    S2_to_1 = Synapses(N2, N1, syn_inter_inh_eqs)
    S2_to_1.connect()
    #S2_to_1.run_regularly()
    S2_to_1.E = Esyn_inh
    S2_to_1.alpha = alpha_inh
    S2_to_1.beta = beta_inh

    # Neuron group state monitors
    M_N1 = StateMonitor(N1, ['x', 'y', 'z'], record=True)
    M_N2 = StateMonitor(N2, 'x', record=True)
    
    run(sim_duration)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    ax1.plot(M_N1.t, M_N1.x[0])
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("x1")
    ax2.plot(M_N2.t, M_N2.x[0])
    ax2.set_xlabel("Time (s)")
    ax2.set_ylabel("x2")
    #plt.savefig("figures/max_noise_60s.png", format="png")
    plt.show()

if __name__ == "__main__":
    main()
