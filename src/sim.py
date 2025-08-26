from brian2 import *
from matplotlib import pyplot as plt
import argparse

def main():
    sim_duration = 120 * second

    # Divide by 1000 to keep unitless without using milli prefix
    Wmax = 0 / 1000000

    tau = 1 * msecond

    pop1_eqs = '''
    dx/dt = (y - a * x ** 3 + b * x ** 2 - z + I_app + C_E * (x_bar - x) +
        sigma * I_syn_intra / amp + sigma * I_syn_inter / amp + xi * Wmax * sqrt(second)) / tau : 1
    dy/dt = (c - d * x ** 2 - y) / tau : 1
    dz/dt = r * (s * (x + x2_bar - x_naught)) / tau : 1    

    a : 1
    b : 1
    c : 1
    d : 1
    s : 1
    r : 1
    x_naught : 1
    I_syn_intra : amp
    I_syn_inter : amp
    I_app : 1
    C_E : 1
    sigma : 1

    x_contrib : 1 (linked)
    x_total = x_contrib : 1 (shared)
    x_bar = x_total / N : 1 (shared)

    # This needs to be passed in throughout runtime
    x2_bar : 1 (shared)
    '''

    pop2_eqs = '''
    dv/dt = (I_app - gL * (v - E_L) - gK * n * (v - E_K) - gCa * m_inf *
        (v - E_Ca) + sigma * C_E * (x_bar - x) * amp + I_syn_inter +
        I_syn_intra - sigma * 0.3 * (z_bar - 3) * amp + sigma * xi * Wmax * sqrt(second)
        * amp) / Cm : volt
    dn/dt = phi * (n_inf - n) / tau_n : 1

    m_inf = 0.5 * (1 + tanh((v - v1) / v2)) : 1
    tau_n = 1 / cosh((v - v3) / (2 * v4)) : 1
    n_inf = 0.5 * (1 + tanh((v - v3) / v4)) : 1

    x = v/(20 * volt) : 1

    Cm : farad
    I_app : amp
    I_syn_inter : amp
    I_syn_intra : amp
    v1 : volt
    v2 : volt
    v3 : volt
    v4 : volt
    phi : Hz
    E_Ca : volt
    E_K : volt
    E_L : volt
    gL : siemens
    gCa : siemens
    gK : siemens
    C_E : 1
    sigma : 1
    z_bar : 1

    x_contrib : 1 (linked)
    x_total = x_contrib : 1 (shared)
    x_bar = x_total / N : 1 (shared)

    x2_bar : 1 (shared)
    '''

    N1 = NeuronGroup(10, pop1_eqs, method='euler')

    N1.x_contrib = linked_var(N1, 'x')

    N1.a = 1
    N1.b = 3
    N1.c = 1
    N1.d = 5
    N1.s = 8
    N1.I_app = 3.1
    N1.x_naught = -4.5
    N1.r = 0.000004
    N1.C_E = 1
    N1.sigma = 1/50

    N2 = NeuronGroup(10, pop2_eqs, method='euler')

    N2.x_contrib = linked_var(N2, 'x')

    N2.Cm = 20 * ufarad
    N2.I_app = 400 * pamp
    N2.v1 = 10 * mvolt
    N2.v2 = 15 * mvolt
    N2.v3 = -1 * mvolt
    N2.v4 = 14.5 * mvolt
    N2.phi = 0.067 * Hz
    N2.E_Ca = 120 * mvolt
    N2.E_K = -84 * mvolt
    N2.E_L = -60 * mvolt
    N2.gL = 2 * msiemens
    N2.gCa = 4 * msiemens
    N2.gK = 8 * msiemens

    universal_syn_eqs ='''
    du/dt = (alpha * T * (1 - u) - beta * u) : 1 (clock-driven)
    T = Tmax / (1 + exp(-(x_pre * mvolt - Vt) / Kp)) : mM

    alpha : mmolar ** -1 * second ** -1
    beta : second ** -1
    '''

    syn_intra_eqs ='''
    I_syn_intra_post = (-G * u * (x_pre * mvolt - E)) : amp (summed)

    G : siemens
    E : volt
    ''' + universal_syn_eqs

    syn_inter_eqs ='''
    I_syn_inter_post = (-G * u * (x_pre * mvolt - E)) : amp (summed)

    #x2_bar_post = x_bar_pre : 1 (summed)
    
    G : siemens
    E : volt
    ''' + universal_syn_eqs

    # Shared synapse parameters
    Tmax = 1 * mM
    Vt = 2 * mV
    Kp = 5 * mV

    # Synapse coupling parameters
    # 1 -> 1, 2 -> 2
    G_intra = [0, 0] * usiemens
    # 1 -> 2, 2 -> 1
    G_inter = [0, 0] * usiemens

    # Excitatory synapse parameters
    alpha_exc = 1.1 / mM / msecond
    beta_exc = 0.19 / msecond
    E_exc = 0 * mvolt

    # Intrapopulation excitatory synapses
    S_exc_intra = Synapses(N1, N1, syn_intra_eqs, method='rk4')
    S_exc_intra.connect(condition='i!=j')
    S_exc_intra.alpha = alpha_exc
    S_exc_intra.beta = beta_exc
    S_exc_intra.E = E_exc
    S_exc_intra.G = G_intra[0]

    # Interpopulation excitatory synapses
    S_exc_inter = Synapses(N1, N2, syn_inter_eqs, method='rk4')
    S_exc_inter.connect(condition='i!=j')
    S_exc_inter.alpha = alpha_exc
    S_exc_inter.beta = beta_exc
    S_exc_inter.E = E_exc
    S_exc_inter.G = G_inter[0]

    # Inhibitory synapse parameters
    alpha_inh = 5 / mM / msecond
    beta_inh = 0.18 / msecond
    E_inh = -80 * mvolt

    # Intrapopulation excitatory synapses
    S_inh_intra = Synapses(N2, N2, syn_intra_eqs, method='rk4')
    S_inh_intra.connect(condition='i!=j')
    S_inh_intra.alpha = alpha_inh
    S_inh_intra.beta = beta_inh
    S_inh_intra.E = E_inh
    S_inh_intra.G = G_intra[1]

    # Interpopulation excitatory synapses
    S_inh_inter = Synapses(N2, N1, syn_inter_eqs, method='rk4')
    S_inh_inter.connect(condition='i!=j')
    S_inh_inter.alpha = alpha_inh
    S_inh_inter.beta = beta_inh
    S_inh_inter.E = E_inh
    S_inh_inter.G = G_inter[1]

    # State monitors
    M_N1 = StateMonitor(N1, ['x', 'y', 'z', 'x_bar'], record=True)
    M_N2 = StateMonitor(N2, ['v', 'x', 'x_bar'], record=True)

    run(sim_duration)

    plt.plot(M_N1.t, M_N2.v[0], label="Population 2")
    #print(M_N1.x[0][0:10])
    plt.xlabel("Time (s)")
    plt.ylabel("v")
    plt.show()

if __name__ == "__main__":
    main()
