"""Plots for the Ca-dependent plasticity gate Wpre on HR->HR synapses.

Recomputes plasticity and sigma_Ca from saved Ca traces so thresholds can be
tweaked without re-running the simulation.
"""
from brian2 import *
import matplotlib.pyplot as plt
import os
import numpy as np
from typing import Dict, Any

import plotting.style  # noqa: F401
from plotting.plotting_utils import *


def plot_plasticity(filepaths: Any, fname: str, params_dict: Dict, title: str,
                    t, x, wpre, u, Ca, plasticity, current, presyn_pop : bool, excite = None, coupling = None,
                    bin_size = 1000) -> None:
    """Five panels: x, Wpre, u, Ca, and plasticity

    Args:
        filepaths (Any): FilePaths with figures_dir.
        fname (str): File name.
        params_dict (Dict): Needs A_LTD, A_LTP, THETA_LTD_START, THETA_LTD_END,
            THETA_LTP_START, CA_SIGMOID_SHIFT, CA_SIGMOID_SLOPE.
        title (str): Title of the plot.
        t: Time vector.
        x: (num_neurons, time) presynaptic population LFP
        wpre: (num_neurons, time) Wpre trace.
        u: (num_neurons, time) u (fraction of open channels) trace.
        Ca: (num_neurons, time) Ca trace.
        plasticity: (num_neurons, time) plasticity trace
        current: (num_neurons, time) synaptic current
        presyn_pop (bool): false is excitatory, true is inhibitory
        binning: bin size in multiples of dt
    """
    has_excite_coupling = int(excite is not None and coupling is not None)
    fig, ax_list = plt.subplots(6 + has_excite_coupling, 1, figsize=(10, 12), sharex=True)
    if has_excite_coupling:
        (ax1, ax2, ax3, ax4, ax5, ax6, ax7) = ax_list
    else:
        (ax1, ax2, ax3, ax4, ax5, ax6) = ax_list
    ax_count = len(ax_list)
    Ca = np.asarray(Ca)
    cb_preval = params_dict['INH_CB_PREVALENCE'] if presyn_pop else params_dict['EXC_CB_PREVALENCE']
    #x_post = np.asarray(x[1])

    effCa = cb_preval * (1 - params_dict['ALPHA_W'] * params_dict['CBD_AMOUNT']) * Ca
    #sigma_Ca = 1 / (1 + exp(-(x_post + params_dict['CA_SIGMOID_SHIFT']) / params_dict['CA_SIGMOID_SLOPE']))

    fig.suptitle(title)
    # Presynaptic LFP
    ax1.plot(t, bin_data(np.mean(x, axis=0), bin_size))
    ax1.set_ylabel("Presynaptic\nLFP (A.U.)", fontsize=10)
    ax1.tick_params(axis='both', which='major', labelsize='10')

    # Synaptic weight
    ax2.plot(t, wpre)
    ax2.set_ylabel("Wpre", fontsize=10)
    ax2.ticklabel_format(useOffset=False)
    ax2.tick_params(axis='both', which='major', labelsize='10')

    # Synaptic opening factor
    ax3.plot(t, bin_data(u, bin_size))
    ax3.set_ylabel("u", fontsize=10)
    ax3.tick_params(axis='both', which='major', labelsize='10')

    # Effective calcium concentration (with dashed lines at plasticity thresholds)
    ax4.plot(t, effCa)
    ax4.set_ylabel("effective Ca", fontsize=10)
    ax4.ticklabel_format(useOffset=False)
    ax4.tick_params(axis='both', which='major', labelsize='10')
    ltd_start_line = ax4.axhline(params_dict['THETA_LTD_START'], 0, 1, c=(0.9, 0.9, 0), dashes=(3,3))
    ltd_end_line = ax4.axhline(params_dict['THETA_LTD_END'], 0, 1, c=(0, 0.9, 0), dashes=(3,3))
    ltp_start_line = ax4.axhline(params_dict['THETA_LTP_START'], 0, 1, c=(0, 0.9, 0.9), dashes=(3,3))
    ax4.legend([ltd_start_line, ltd_end_line, ltp_start_line], ['LTD start', 'LTD end', 'LTP start'], loc='upper right', fontsize=8)

    # Plasticity
    ax5.plot(t, bin_data(plasticity, bin_size))
    ax5.set_ylabel("Plasticity", fontsize=10)
    ax5.tick_params(axis='both', which='major', labelsize=10)

    # Synaptic current
    ax6.plot(t, bin_data(np.mean(current, axis=0), bin_size))
    ax6.set_ylabel("Synaptic Current", fontsize=10)
    ax6.yaxis.get_offset_text().set_fontsize(10)
    ax6.tick_params(axis='both', which='major', labelsize=10)

    # Coupling and Epop excitability
    if has_excite_coupling:
        ax7.plot(t, excite, label='x0', color='blue')
        ax7b = ax7.twinx()
        ax7b.plot(t, coupling, label='Ce', color='orange')
        ax7.set_ylabel("Epop excite", fontsize=10, c='blue')
        ax7b.set_ylabel("Coupling (J)", fontsize=10, c='orange')
        ax7.tick_params(axis='y', labelcolor='blue', labelsize=10)
        ax7b.tick_params(axis='y', labelcolor='orange', labelsize=10)

    # Intra- and interpopulation synapse strength

    ax_list[ax_count - 1].set_xlabel("time (s)", fontsize=10)
    ax_list[ax_count - 1].tick_params(axis='x', labelsize=10)
    fig.tight_layout()

    plt.savefig(os.path.join(filepaths.figures_dir, fname), format="png")
    plt.show()
