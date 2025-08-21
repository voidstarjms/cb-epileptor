from brian2 import *
from matplotlib import pyplot as plt
import argparse

def main():
    pop1_eqs = '''
    dx/dt = y - a * x ** 3 + b * x ** 2 - z + I_app + C_E * (x_bar - x) +
        sigma * I_syn_intra + sigma * I_syn_inter + W(t) : 1
    dy/dt = c - d * x ** 2 - y : 1
    dz/dt = r * (s * (x + x2_bar - x_naught)) : 1

    a : 1
    b : 1
    c : 1
    d : 1
    s : 1
    I : 1
    r : 1
    x_naught : 1
    I_syn_intra : 1
    I_syn_inter : 1
    I_app : 1
    C_E : 1
    sigma : 1

    # This needs to be passed in throughout runtime
    x2_bar : 1
    '''

    pop2_eqs = '''
    dv/dt = (I_app - gL * (v - E_L) - gK * n * (V - E_K) - gCa * m_inf *
        (v - E_Ca) + sigma * C_E * (x_bar - x) + I_syn_inter + I_syn_intra -
        sigma * 0.3 * (z_bar - 3) + sigma * W(t)) / Cm : volt
    dn/dt = phi * (n_inf - n) / tau_n : 1

    m_inf = 0.5 * (1 + tanh((v - v1) / v2)) : 1
    tau_n = 1 / cosh((v - v3) / (2 * v4)) : 1
    n_inf = 0.5 * (1 + tanh((v - v3) / v4)) : 1

    x = v/20 : 1

    Cm : farad/meter ** 2
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
    C_E : 1
    sigma : 1
    '''

    N1 = NeuronGroup(10, pop1_eqs, method='rk4')

    N1.x = 0
    N1.y = 0
    N1.z = 0
    N1.a = 1
    N1.b = 3
    N1.c = 1
    N1.d = 5
    N1.s = 8
    N1.I_app = 3.1
    N1.x_naught = -4.5
    N1.r = 0.000004
    N1.C_E = 0
    N1.sigma = 1/50

    N2 = NeuronGroup(10, pop2_eqs, method='rk4')

    N2.Cm = 20 * ufarad / cmeter ** 2
    N2.I_app = 40 * amp
    N2.v1 = 1.2 * volt
    N2.v2 = 18 * volt
    N2.v3 = 12 * volt
    N2.v4 = 17.4 * volt
    N2.phi = 0.067 * Hz
    N2.E_Ca = 120 * volt
    N2.E_K = -84 * volt
    N2.E_L = -60 * volt

    universal_syn_eqs ='''
    du/dt = alpha * T * (1 - u) - beta * u : 1 (clock-driven)
    T = Tmax / (1 + exp(-(x_pre - Vt) / Kp)) : 1

    alpha : 1
    beta : 1
    Tmax : 1
    Vt : 1
    Kp : 1
    '''

    syn_intra_eqs ='''
    I_syn_intra_post = -G * u * (x_pre - E) : 1

    G : 1
    E : 1
    ''' + universal_syn_eqs

    syn_inter_eqs ='''
    I_syn_inter_post = -G * u * (x_pre - E) : 1

    G : 1
    E : 1
    ''' + universal_syn_eqs

    # Shared synapse parameters
    Tmax = 1
    Vt = 2
    Kp = 5

    # Synapse coupling parameters
    G_intra = [0, 0.4]
    G_inter = [0, 0.4]

    # Excitatory synapse parameters
    alpha_exc = 1.1
    beta_exc = 0.19
    E_exc = 0

    # Intrapopulation excitatory synapses
    S_exc_intra = Synapses(N1, N1, syn_intra_eqs, method='rk4')
    S_exc_intra.connect(condition='i!=j')
    S_exc_intra.alpha = alpha_exc
    S_exc_intra.beta = beta_exc
    S_exc_intra.Tmax = Tmax
    S_exc_intra.E = E_exc
    S_exc_intra.G = G_intra[0]
    S_exc_intra.Vt = Vt
    S_exc_intra.Kp = Kp

    # Interpopulation excitatory synapses
    S_exc_inter = Synapses(N1, N2, syn_inter_eqs, method='rk4')
    S_exc_inter.connect(condition='i!=j')
    S_exc_inter.alpha = alpha_exc
    S_exc_inter.beta = beta_exc
    S_exc_inter.Tmax = Tmax
    S_exc_inter.E = E_exc
    S_exc_inter.G = G_inter[0]
    S_exc_inter.Vt = Vt
    S_exc_inter.Kp = Kp

    # Inhibitory synapse parameters
    alpha_inh = 5
    beta_inh = 0.18
    E_inh = -80

    # Intrapopulation excitatory synapses
    S_inh_intra = Synapses(N2, N2, syn_intra_eqs, method='rk4')
    S_inh_intra.connect(condition='i!=j')
    S_inh_intra.alpha = alpha_inh
    S_inh_intra.beta = beta_inh
    S_inh_intra.Tmax = Tmax
    S_inh_intra.E = E_inh
    S_inh_intra.G = G_intra[1]
    S_inh_intra.Vt = Vt
    S_inh_intra.Kp = Kp

    # Interpopulation excitatory synapses
    S_inh_inter = Synapses(N2, N1, syn_inter_eqs, method='rk4')
    S_inh_inter.connect(condition='i!=j')
    S_inh_inter.alpha = alpha_inh
    S_inh_inter.beta = beta_inh
    S_inh_inter.Tmax = Tmax
    S_inh_inter.E = E_inh
    S_inh_inter.G = G_inter[1]
    S_inh_inter.Vt = Vt
    S_inh_inter.Kp = Kp

    # State monitors
    M_N1 = StateMonitor(N1, ['x', 't'], record=True)
    M_N2 = StateMonitor(N2, 'v', record=True)

    run(10 * second)

    plt.plot(M_N1.t, M_N1.x, label="Population 1")
    plt.show()

if __name__ == "__main__":
    main()
