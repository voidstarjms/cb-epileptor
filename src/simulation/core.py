from brian2 import *
from brian2tools import *
import numpy as np
import params
import data_processing
def run_sim(cb_on=True, data_dir=None, g_intra_override=None, g_inter_override=None):
    # Setup Simulation
    defaultclock.dt = params.TAU_CLOCK / params.DT_SCALING
    print("defaultclock.dt is: ", defaultclock.dt)

    # Build timed arrays from current params at call time
    timed_x_naught = TimedArray(params.X_NAUGHT_VALS, dt=params.SIM_DURATION // len(params.X_NAUGHT_VALS))
    timed_coupling_strength = TimedArray(params.COUPLING_VALS, dt=params.SIM_DURATION // len(params.COUPLING_VALS))
    timed_G_inter = TimedArray(params.G_INTER_VALS, dt=params.SIM_DURATION // len(params.G_INTER_VALS))
    timed_G_intra = TimedArray(params.G_INTRA_VALS, dt=params.SIM_DURATION // len(params.G_INTRA_VALS))

    # --- Population 1: Hindmarsh-Rose ---
    pop1_eqs = '''
    dx/dt = (y - a * x ** 3 + b * x ** 2 - z + I_app_1
        + ISOLATE * (timed_CE(t) * (x_bar - x)
        + Wmax * xi * sqrt(second)
        + sigma_1 * (I_syn_intra + I_syn_inter) / I_scale)) / tau : 1
    dy/dt = (c - d * x ** 2 - y) / tau : 1
    dz/dt = r * (s * (x + ISOLATE * x2_bar - timed_x0(t)) - z_bar) : 1

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
        'ISOLATE': params.ISOLATE,
        'Wmax': params.W_MAX,
        'I_scale': params.I_SCALE,
        'timed_x0': timed_x_naught, 'timed_CE': timed_coupling_strength,
    }

    # CHANGE HOW X0 IS LOADED - SHOULD USE TIMED ARRAY[0]
    N1 = NeuronGroup(params.NUM_CELLS, pop1_eqs, method='euler',
                     threshold=params.HR_THRESHOLD, reset='',
                     namespace=pop1_namespace, refractory=params.HR_REFRACTORY_CONDITION)

    N1.x = np.ones(params.NUM_CELLS) * (params.HR_X_NAUGHT+1.5) + randn(params.NUM_CELLS) * params.W_MAX
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
        + timed_CE(t) * (x_bar - x) - ml_z_bar_scale * (z_bar - ml_z_bar_offset))
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
        'Wmax': params.W_MAX,
        'ISOLATE': params.ISOLATE,
        'timed_CE': timed_coupling_strength,
        'ml_z_bar_scale': params.ML_Z_BAR_SCALE,
        'ml_z_bar_offset': params.ML_Z_BAR_OFFSET,
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
        'Kp': params.SYN_KP,
        'tau_wpre': params.TAU_WPRE,
        'tau_ca': params.TAU_CA,
        'theta_ltd_start': params.THETA_LTD_START,
        'theta_ltd_end': params.THETA_LTD_END,
        'theta_ltp_start': params.THETA_LTP_START,
        'A_ltp': params.A_LTP,
        'A_ltd': params.A_LTD,
        'timed_G_intra': g_intra_override if g_intra_override is not None else timed_G_intra,
        'timed_G_inter': g_inter_override if g_inter_override is not None else timed_G_inter,
        'ca_sigmoid_shift': params.CA_SIGMOID_SHIFT,
        'ca_sigmoid_slope': params.CA_SIGMOID_SLOPE,
    }

    syn_input_scale = 1/pop1_namespace['sigma_1']

    syn_eqs_pre = '''
        du/dt = (alpha * T * (1 - u) - beta * u) : 1 (clock-driven)
        T = Tmax / (1 + exp(-(x_pre * syn_input_scale * mvolt - Vt) / Kp)) : mM

        plasticity = 1 - A_ltd * int(Ca > theta_ltd_start) * int(Ca < theta_ltd_end) + A_ltp * int(Ca > theta_ltp_start) : 1
        dWpre/dt = (plasticity - Wpre) / tau_wpre : 1 (clock-driven)
        dCa/dt = (sigma_Ca - Ca) / tau_ca : 1 (clock-driven)
        sigma_Ca = 1 / (1 + exp(-(x_post + ca_sigmoid_shift) / ca_sigmoid_slope)) : 1

        E : volt
        alpha : mmolar ** -1 * second ** -1
        beta : second ** -1
    '''


    wpre_term = 'Wpre' if cb_on else '1'

    # S2_to_2: pre=Pop2, post=Pop2
    intra_syn_eqs = f'''
    I_syn_intra_post = (-timed_G_intra(t) * u * (x_post * syn_input_scale * mvolt - E)) * {wpre_term} : amp (summed)
    ''' + syn_eqs_pre

    # S1_to_2: pre=Pop1, post=Pop2
    inter_syn_eqs = f'''
    I_syn_inter_post = (-timed_G_inter(t) * u * (x_post * syn_input_scale * mvolt - E)) * {wpre_term} : amp (summed)
    ''' + syn_eqs_pre



    S1_to_1 = Synapses(N1, N1, intra_syn_eqs, method='euler', namespace=syn_namespace)
    S1_to_1.connect()
    S1_to_1.E = params.SYN_E_EXC
    S1_to_1.alpha = params.SYN_ALPHA_EXC
    S1_to_1.beta = params.SYN_BETA_EXC

    S1_to_2 = Synapses(N1, N2, inter_syn_eqs, method='euler', namespace=syn_namespace)
    S1_to_2.connect()
    S1_to_2.run_regularly('z_bar_post = z_bar_pre', dt=defaultclock.dt)
    S1_to_2.E = params.SYN_E_EXC
    S1_to_2.alpha = params.SYN_ALPHA_EXC
    S1_to_2.beta = params.SYN_BETA_EXC

    S2_to_2 = Synapses(N2, N2, intra_syn_eqs, method='euler', namespace=syn_namespace)
    S2_to_2.connect()
    S2_to_2.E = params.SYN_E_INH
    S2_to_2.alpha = params.SYN_ALPHA_INH
    S2_to_2.beta = params.SYN_BETA_INH

    S2_to_1 = Synapses(N2, N1, inter_syn_eqs, method='euler', namespace=syn_namespace)
    S2_to_1.connect()
    S2_to_1.run_regularly('x2_bar_post = x_bar_pre', dt=defaultclock.dt)
    S2_to_1.E = params.SYN_E_INH
    S2_to_1.alpha = params.SYN_ALPHA_INH
    S2_to_1.beta = params.SYN_BETA_INH

    # Don't record transient period
    run(params.TRANSIENT*second)

    M_N1 = StateMonitor(N1, ['x', 'y', 'z', 'I_syn_inter', 'I_syn_intra'], record=True)
    M_N2 = StateMonitor(N2, ['x', 'n', 'I_syn_inter'], record=True)

    SM_N1 = SpikeMonitor(N1)
    SM_N2 = SpikeMonitor(N2)

    M_S1_1 = StateMonitor(S1_to_1, ['Wpre', 'Ca', 'u'], record=True)

    # Run
    run(params.SIM_DURATION)
    data_processing.save_data(M_N1, M_N2, SM_N1, SM_N2, M_S1_1, cb_on, data_dir=data_dir)
