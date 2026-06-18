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

def plot_plasticity(filepaths: Any, fname: str, params_dict: Dict, title: str, t, x, wpre, u, Ca) -> None:
    """Five panels: x, Wpre, u, Ca, and plasticity

    Plasticity and sigma_Ca are recomputed here so thresholds can be tweaked
    without re-running the sim.

    Args:
        filepaths (Any): FilePaths with figures_dir.
        fname (str): File name.
        params_dict (Dict): Needs A_LTD, A_LTP, THETA_LTD_START, THETA_LTD_END,
            THETA_LTP_START, CA_SIGMOID_SHIFT, CA_SIGMOID_SLOPE.
        title (str): Title of the plot.
        t: Time vector.
        x: (num_neurons, time) x
        wpre: (num_neurons, time) Wpre trace.
        u: (num_neurons, time) u (fraction of open channels) trace.
        Ca: (num_neurons, time) Ca trace.
    """
    fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, figsize=(10, 8), sharex=True)
    Ca = np.asarray(Ca)
    #x_post = np.asarray(x[1])

    plasticity = (1
        - params_dict['A_LTD'] * np.where(Ca > params_dict['THETA_LTD_START'], 1, 0)
        * np.where(Ca < params_dict['THETA_LTD_END'], 1, 0)
        + params_dict['A_LTP'] * (Ca > params_dict['THETA_LTP_START']))
    #sigma_Ca = 1 / (1 + exp(-(x_post + params_dict['CA_SIGMOID_SHIFT']) / params_dict['CA_SIGMOID_SLOPE']))


    fig.suptitle(title)
    ax1.plot(t, np.mean(x))
    ax1.set_ylabel("Presynaptic LFP (A.U.)")
    ax2.plot(t, np.mean(wpre))
    ax2.set_ylabel("Wpre")
    ax3.plot(t, np.mean(u))
    ax3.set_ylabel("u")
    ax4.plot(t, np.mean(Ca))
    ax4.set_ylabel("Ca")
    ax5.plot(t, np.mean(plasticity))
    ax5.set_ylabel("Plasticity")
    #ax6.plot(t, sigma_Ca)
    #ax6.set_ylabel("ca_signal")
    ax5.set_xlabel("time (s)")

    plt.savefig(os.path.join(filepaths.figures_dir, fname), format="png")
    plt.show()
