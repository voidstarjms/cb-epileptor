"""Brian2 simulation of the two-population seizure model (Naze et al. 2015).

Builds two coupled populations (Hindmarsh-Rose excitatory + Morris-Lecar
inhibitory), wires up intra/inter-population chemical synapses with calcium
control plasticity, and runs the simulation. Saves output via data_processing.
"""
from brian2 import *
import numpy as np
import data_processing
from typing import Dict, Any

def run_sim(filepaths: Any, params_dict: Dict[str, Any], cb_on: bool = True, save: bool = True) -> None:
    """Build both populations and their synapses, run the sim, save the output.

    If cb_on is False, Wpre is dropped from the synaptic current (plasticity
    still evolves but has no effect on the dynamics).

    Args:
        filepaths (Any): FilePaths. save_data writes to filepaths.data_dir/output.pkl.
        params_dict (Dict[str, Any]): Flat dict from param_loader. Expected keys
            listed in run.REQUIRED_PARAMS.
        cb_on (bool): If False, Wpre is held at 1 in the synaptic current.
    """
    # Setup Simulation
    defaultclock.dt = params_dict['TAU_CLOCK'] / params_dict['DT_SCALING']
    print("defaultclock.dt is: ", defaultclock.dt)

    print("============================================In model=============================================")
    print(params_dict)

    # Simulation-control namespace hyper-params (used throughout)
    sim_namespace = {
        'num_cells':    params_dict['NUM_CELLS'],
        'sim_duration': params_dict['SIM_DURATION'],
        'transient':    params_dict['TRANSIENT'],
    }

    # Time-varying schedules. Each VALS list is stepped through evenly across SIM_DURATION.
    timed_x_naught = TimedArray(params_dict['BASE_EXCITE_VALS'], dt=params_dict['BASE_EXCITE_DT'])
    timed_coupling_strength = TimedArray(params_dict['COUPLING_VALS'], dt=params_dict['COUPLING_DT'])

    # compute effective synapse strength when factoring network connection percentage
    eff_g_inter_vals = params_dict['G_INTER_VALS'] / params_dict['PCT_CONNECT']
    eff_g_intra_vals = params_dict['G_INTRA_VALS'] / params_dict['PCT_CONNECT']

    timed_G_inter = TimedArray(eff_g_inter_vals, dt=params_dict['G_INTER_DT'])
    timed_G_intra = TimedArray(eff_g_intra_vals, dt=params_dict['G_INTRA_DT'])
    
    # --- Population 1: Hindmarsh-Rose ---
    pop1_namespace = {
        'a': params_dict['EPOP_A'], 'b': params_dict['EPOP_B'], 'c': params_dict['EPOP_C'],
        'd': params_dict['EPOP_D'], 's': params_dict['EPOP_S'], 'I_app_1': params_dict['EPOP_I_APP'],
        'x_naught': params_dict['EPOP_BASE_EXCITE'], 'r': params_dict['EPOP_R'],
        'sigma_1': params_dict['EPOP_SIGMA'], 'tau': params_dict['TAU_CLOCK'],
        'ISOLATE': params_dict['ISOLATE'], 'NOISE_INIT_OFFSET': params_dict['NOISE_INIT_OFFSET'],
        'Wmax': params_dict['W_MAX'],
        'I_scale': params_dict['I_SCALE'],
        'J': params_dict['COUPLING_STRENGTH'],
    }

    pop1_eqs = '''
    dx/dt = (y - a * x ** 3 + b * x ** 2 - z + I_app_1
        + ISOLATE * (timed_coupling_strength(t) * (x_bar - x)
        + Wmax * xi * sqrt(second)
        + sigma_1 * (I_syn_intra + I_syn_inter) / I_scale)) / tau : 1
    dy/dt = (c - d * x ** 2 - y) / tau : 1
    dz/dt = r * (s * (x + ISOLATE * x2_bar - timed_x_naught(t)) - z_bar) : 1

    x_bar : 1
    z_bar : 1
    x2_bar : 1
    I_syn_intra : amp
    I_syn_inter : amp
    x0_t = timed_x_naught(t) : 1
    ce_t = timed_coupling_strength(t) : 1
    '''

    N1 = NeuronGroup(sim_namespace['num_cells'], pop1_eqs, method='euler',
                     threshold=params_dict['EPOP_THRESHOLD'], reset='',
                     namespace=pop1_namespace, refractory=params_dict['EPOP_REFRACTORY_CONDITION'])

    N1.x = np.ones(sim_namespace['num_cells']) * (timed_x_naught(0 * second)+ pop1_namespace['NOISE_INIT_OFFSET']) + randn(sim_namespace['num_cells']) * pop1_namespace['Wmax']
    N1.y = 'c - d*x**2'
    N1.z = '(s*(x - timed_x_naught(0 * second)))'

    # Population 1 Averaging
    exc_averaging_eqs ='''
    x_bar_post = x_pre / num_cells : 1 (summed)
    z_bar_post = z_pre / num_cells : 1 (summed)
    '''
    gap_junctions_1 = Synapses(N1, N1, exc_averaging_eqs, namespace={'num_cells': sim_namespace['num_cells']})
    gap_junctions_1.connect()


    # --- Population 2: Morris-Lecar ---
    pop2_namespace = {
        'Cm': params_dict['IPOP_CM'], 'I_app_2': params_dict['IPOP_I_APP'], 'gL': params_dict['IPOP_GL'],
        'E_L': params_dict['IPOP_E_L'], 'gK': params_dict['IPOP_GK'], 'E_K': params_dict['IPOP_E_K'],
        'gCa': params_dict['IPOP_GCA'], 'E_Ca': params_dict['IPOP_E_CA'],
        'h_Ca': params_dict['IPOP_H_CA'], 'lambda_Ca': params_dict['IPOP_LAMBDA_CA'], 'h_K': params_dict['IPOP_H_K'], 'lambda_K': params_dict['IPOP_LAMBDA_K'],
        'phi': params_dict['IPOP_PHI'], 'sigma_2': params_dict['IPOP_SIGMA'],
        'Wmax': params_dict['W_MAX'],
        'ISOLATE': params_dict['ISOLATE'],
        'J': params_dict['COUPLING_STRENGTH'], 
        'ml_z_bar_scale': params_dict['IPOP_Z_BAR_SCALE'],
        'ml_z_bar_offset': params_dict['IPOP_Z_BAR_OFFSET'],
    }

    pop2_eqs = '''
    dv/dt = (I_app_2 - gL*(v-E_L) - gK*n*(v-E_K) - gCa*m_inf*(v-E_Ca)
        + ISOLATE * (sigma_2 * (Wmax * xi * sqrt(second)
        + timed_coupling_strength(t) * (x_bar - x) - ml_z_bar_scale * (z_bar - ml_z_bar_offset))
        + (I_syn_intra + I_syn_inter))) / Cm : volt
    dn/dt = phi * (n_inf - n) / tau_n : 1

    m_inf = 0.5 * (1 + tanh((v - h_Ca) / lambda_Ca)) : 1
    n_inf = 0.5 * (1 + tanh((v - h_K) / lambda_K)) : 1
    tau_n = 1 / cosh((v - h_K) / (2 * lambda_K)) : 1

    x = v/(20 * mV) : 1

    x_bar : 1
    z_bar : 1
    I_syn_inter : amp
    I_syn_intra : amp
    '''

    N2 = NeuronGroup(sim_namespace['num_cells'], pop2_eqs, method='euler',
                     threshold=params_dict['IPOP_THRESHOLD'], reset='',
                     namespace=pop2_namespace, refractory=params_dict['IPOP_REFRACTORY_CONDITION'])

    N2.v = pop2_namespace['E_L'] * np.ones(sim_namespace['num_cells']) + \
           randn(sim_namespace['num_cells']) * pop2_namespace['Wmax'] * volt
    N2.n = 'n_inf'

    # Population 2 Averaging
    inh_averaging_eqs ='''
    x_bar_post = x_pre / num_cells : 1 (summed)
    '''
    gap_junctions_2 = Synapses(N2, N2, inh_averaging_eqs, namespace={'num_cells': sim_namespace['num_cells']})
    gap_junctions_2.connect()


    # --- Synapses ---
    syn_namespace = {
        'pct_connect': params_dict['PCT_CONNECT'],
        'Tmax': params_dict['SYN_TMAX'],
        'Vt': params_dict['SYN_VT'],
        'Kp': params_dict['SYN_KP'],
        'tau_wpre': params_dict['TAU_WPRE'],
        'tau_ca': params_dict['TAU_CA'],
        'theta_ltd_start': params_dict['THETA_LTD_START'],
        'theta_ltd_end': params_dict['THETA_LTD_END'],
        'theta_ltp_start': params_dict['THETA_LTP_START'],
        'A_ltp': params_dict['A_LTP'],
        'A_ltd': params_dict['A_LTD'],
        'timed_G_intra': timed_G_intra,
        'timed_G_inter': timed_G_inter,
        'ca_sigmoid_shift': params_dict['CA_SIGMOID_SHIFT'],
        'ca_sigmoid_slope': params_dict['CA_SIGMOID_SLOPE'],
        'E_exc': params_dict['SYN_E_EXC'],
        'alpha_exc': params_dict['SYN_ALPHA_EXC'],
        'beta_exc': params_dict['SYN_BETA_EXC'],
        'E_inh': params_dict['SYN_E_INH'],
        'alpha_inh': params_dict['SYN_ALPHA_INH'],
        'beta_inh': params_dict['SYN_BETA_INH'],
        'alpha_w': params_dict['ALPHA_W'],
        'beta_ca': params_dict['BETA_CA'],
        'gamma_ca': params_dict['GAMMA_CA'],
        'cbd': params_dict['CBD_AMOUNT'],
    }

    syn_input_scale = 1/pop1_namespace['sigma_1']

    syn_eqs_pre = '''
        du/dt = (alpha * T * (1 - u) - beta * u) : 1 (clock-driven)
        T = Tmax / (1 + exp(-(x_pre * syn_input_scale * mvolt - Vt) / Kp)) : mM

        effCa = (1 - alpha_w * cbd) * Ca : 1
        plasticity = 1 - A_ltd * int(effCa > theta_ltd_start) * int(effCa < theta_ltd_end) + A_ltp * int(effCa > theta_ltp_start) : 1
        dWpre/dt = (plasticity - Wpre) / tau_wpre : 1 (clock-driven)
        dCa/dt = (sigma_Ca - (beta_ca - gamma_ca * cbd) * Ca) / tau_ca : 1 (clock-driven)
        sigma_Ca = 1 / (1 + exp(-(x_post + ca_sigmoid_shift) / ca_sigmoid_slope)) : 1

        E : volt
        alpha : mmolar ** -1 * second ** -1
        beta : second ** -1
    '''

    wpre_term = 'Wpre' if cb_on else '1'

    # S1_to_1: pre=Pop1, post=Pop1
    intra_syn_eqs = f'''
    I_syn_intra_post = (-timed_G_intra(t) * u * (x_post * syn_input_scale * mvolt - E)) * {wpre_term} : amp (summed)
    ''' + syn_eqs_pre

    # S1_to_2: pre=Pop1, post=Pop2
    inter_syn_eqs = f'''
    I_syn_inter_post = (-timed_G_inter(t) * u * (x_post * syn_input_scale * mvolt - E)) * {wpre_term} : amp (summed)
    ''' + syn_eqs_pre

    # Exc-to-Exc synapses
    S1_to_1 = Synapses(N1, N1, intra_syn_eqs, method='euler', namespace=syn_namespace)
    S1_to_1.connect(p=syn_namespace['pct_connect'])
    S1_to_1.E = syn_namespace['E_exc']
    S1_to_1.alpha = syn_namespace['alpha_exc']
    S1_to_1.beta = syn_namespace['beta_exc']

    # Exc-to-Inh synapses
    S1_to_2 = Synapses(N1, N2, inter_syn_eqs, method='euler', namespace=syn_namespace)
    S1_to_2.connect(p=syn_namespace['pct_connect'])
    # Pass HR's z_bar into ML so pop2's eqs can read it.
    S1_to_2.run_regularly('z_bar_post = z_bar_pre', dt=defaultclock.dt)
    S1_to_2.E = syn_namespace['E_exc']
    S1_to_2.alpha = syn_namespace['alpha_exc']
    S1_to_2.beta = syn_namespace['beta_exc']

    # Inh-to-Inh synapses
    S2_to_2 = Synapses(N2, N2, intra_syn_eqs, method='euler', namespace=syn_namespace)
    S2_to_2.connect(p=syn_namespace['pct_connect'])
    S2_to_2.E = syn_namespace['E_inh']
    S2_to_2.alpha = syn_namespace['alpha_inh']
    S2_to_2.beta = syn_namespace['beta_inh']

    # Inh-to-Exc
    S2_to_1 = Synapses(N2, N1, inter_syn_eqs, method='euler', namespace=syn_namespace)
    S2_to_1.connect(p=syn_namespace['pct_connect'])
    # Pass ML's x_bar into HR as x2_bar so pop1's dz/dt can read it.
    S2_to_1.run_regularly('x2_bar_post = x_bar_pre', dt=defaultclock.dt)
    S2_to_1.E = syn_namespace['E_inh']
    S2_to_1.alpha = syn_namespace['alpha_inh']
    S2_to_1.beta = syn_namespace['beta_inh']

    # Run the transient before attaching monitors so it isn't recorded.
    run(sim_namespace['transient']*second)
    print(sim_namespace['transient'])

    # Population state monitors
    M_N1 = StateMonitor(N1, ['x', 'y', 'z', 'I_syn_inter', 'I_syn_intra'], record=True)
    M_N2 = StateMonitor(N2, ['x', 'n', 'I_syn_inter'], record=True)
    # Parameter state monitor
    M_PARAM = StateMonitor(N1, ['x0_t', 'ce_t'], record=0)
    # Population spike monitors
    SM_N1 = SpikeMonitor(N1)
    SM_N2 = SpikeMonitor(N2)
    # Synaptic plasticity state monitors
    M_S1_1 = StateMonitor(S1_to_1, ['Wpre', 'Ca', 'u'], record=True)
    M_S1_2 = StateMonitor(S1_to_2, ['Wpre', 'Ca', 'u'], record=True)
    M_S2_2 = StateMonitor(S2_to_2, ['Wpre', 'Ca', 'u'], record=True)
    M_S2_1 = StateMonitor(S2_to_1, ['Wpre', 'Ca', 'u'], record=True)
    
    # Run
    run(sim_namespace['sim_duration'])
    sim_data = data_processing.pack_data(params_dict, M_N1, M_N2, SM_N1, SM_N2, M_S1_1, M_S1_2, M_S2_1, M_S2_2, cb_on=cb_on, M_param=M_PARAM)
    if save == True:
        data_processing.save_data(filepaths, sim_data)
    return sim_data
