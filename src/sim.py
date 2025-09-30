from brian2 import *
from matplotlib import pyplot as plt
import argparse

def main():
    sim_duration = 60 * second

    tau = 1 * msecond

    num_cells = 1

    Wmax = 0

    ### Population 1 (Hindmarsh-Rose)

    # Population 1 parameters
    a = 1
    b = 3
    c = 1
    d = 5
    s = 8
    I_app_1 = 3.1
    x_naught = -1.6
    r = 0.000004

    # Population 1 equations
    pop1_eqs = '''
    dx/dt = (y - a * x ** 3 + b * x ** 2 - z + I_app_1 + sigma_1 * (x_bar - x)
        + sigma_1 * (2 * Wmax * xi * sqrt(second) - Wmax)) / tau : 1
    dy/dt = (c - d * x ** 2 - y) / tau : 1
    dz/dt = r * (s * (x - x_naught) - z) / tau : 1

    x_bar : 1    
    '''

    N1 = NeuronGroup(num_cells, pop1_eqs, method='euler') 
    N1.x = x_naught
    N1.y = 'c - d*x**2'
    N1.z = 'r*(s*(x - x_naught))'

    # Population 1 gap junctions
    sigma_1 = 1/50
    gap_junction_eqs ='''
    x_bar_post = x / num_cells : 1 (summed)
    '''
    gap_junctions_1 = Synapses(N1, N1, gap_junction_eqs)

    ### Population 2 (Morris-Lecar)

    # Population 2 parameters
    Cm = 20 * ufarad
    I_app_2 = 400 * uamp
    v1 = 1.2 * mvolt
    v2 = 18 * mvolt
    v3 = 12 * mvolt
    v4 = 17.4 * mvolt
    phi = 0.067 * Hz
    E_Ca = 120 * mvolt
    E_K = -84 * mvolt
    E_L = -60 * mvolt
    gL = 2 * msiemens
    gCa = 4 * msiemens
    gK = 8 * msiemens

    # Population 2 equations    
    pop2_eqs = '''
    dv/dt = (I_app_2 - gL * (v - E_L) - gK * n * (v - E_K) - gCa * m_inf *
        (v - E_Ca) + sigma_2 * (2 * Wmax * xi * sqrt(second) - Wmax) +
        sigma_2 * (x_bar - x)) / Cm : volt
    dn/dt = phi * (n_inf - n) / tau_n : 1

    m_inf = 0.5 * (1 + tanh((v - v1) / v2)) : 1
    tau_n = 1 / cosh((v - v3) / (2 * v4)) : 1
    n_inf = 0.5 * (1 + tanh((v - v3) / v4)) : 1

    x = v/(20 * mV) : 1
    
    x_bar : 1
    '''
    
    N2 = NeuronGroup(num_cells, pop2_eqs, method='euler')
    N2.v = E_L
    N2.n = 'n_inf'
    
    # Population 2 gap junctions
    sigma_2 = 50 * uA
    gap_junctions_2 = Synapses(N2, N2, gap_junction_eqs)

    # universal_syn_eqs ='''
    # du/dt = (alpha * T * (1 - u) - beta * u) : 1 (clock-driven)
    # T = Tmax / (1 + exp(-(x_pre * mvolt - Vt) / Kp)) : mM

    # alpha : mmolar ** -1 * second ** -1
    # beta : second ** -1
    # '''

    # syn_intra_eqs ='''
    # I_syn_intra_post = (-G * u * (x_post * mvolt - E)) : amp (summed)

    # G : siemens
    # E : volt
    # ''' + universal_syn_eqs

    # syn_inter_eqs ='''
    # I_syn_inter_post = (-G * u * (x_post * mvolt - E)) : amp (summed)

    # G : siemens
    # E : volt
    # ''' + universal_syn_eqs

    # Neuron group state monitors
    M_N1 = StateMonitor(N1, ['x', 'y', 'z'], record=True)
    M_N2 = StateMonitor(N2, 'x', record=True)
    
    run(sim_duration)

    print(M_N1.x[0][0:10])
    print(M_N1.y[0][0:10])
    print(M_N1.z[0][0:10])

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
