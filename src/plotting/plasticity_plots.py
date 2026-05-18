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

def plot_wpre(filepaths: Any, params_dict: Dict, t, x, wpre, u, Ca) -> None:
    """Seven panels: x of two neurons, Wpre, u, Ca, plasticity, and sigma_Ca.

    Plasticity and sigma_Ca are recomputed here so thresholds can be tweaked
    without re-running the sim.

    Args:
        filepaths (Any): FilePaths with figures_dir.
        params_dict (Dict): Needs A_LTD, A_LTP, THETA_LTD_START, THETA_LTD_END,
            THETA_LTP_START, CA_SIGMOID_SHIFT, CA_SIGMOID_SLOPE.
        t: Time vector.
        x: (num_neurons, time) HR x. Uses x[0] and x[1].
        wpre: Wpre trace for one synapse.
        u: u (neurotransmitter) trace for one synapse.
        Ca: Ca trace for one synapse.
    """
    fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7) = plt.subplots(7, 1, figsize=(10, 8), sharex=True)
    Ca = np.asarray(Ca)
    x_post = np.asarray(x[1])

    plasticity = (1
        - params_dict['A_LTD'] * (Ca > params_dict['THETA_LTD_START']).astype(float) * (Ca < params_dict['THETA_LTD_END']).astype(float)
        + params_dict['A_LTP'] * (Ca > params_dict['THETA_LTP_START']).astype(float))
    sigma_Ca = 1 / (1 + exp(-(x_post + params_dict['CA_SIGMOID_SHIFT']) / params_dict['CA_SIGMOID_SLOPE']))


    fig.suptitle("Wpre Within Excitatory Population")
    ax1.plot(t, x[0])
    ax1.set_ylabel("Neuron 1 x ")
    ax2.plot(t, x[1])
    ax2.set_ylabel("Neuron 2 x ")
    ax3.plot(t, wpre)
    ax3.set_ylabel("Wpre")
    ax4.plot(t, u)
    ax4.set_ylabel("u")
    ax5.plot(t, Ca)
    ax5.set_ylabel("Ca")
    ax6.plot(t, plasticity)
    ax6.set_ylabel("plasticity")
    ax7.plot(t, sigma_Ca)
    ax7.set_ylabel("ca_signal")
    ax7.set_xlabel("time (s)")

    plt.savefig(os.path.join(filepaths.figures_dir, "N1_to_1_wpre.png"), format="png")
    plt.show()
