from brian2 import *
import matplotlib.pyplot as plt
import os
import numpy as np
from typing import Dict, Any

def plot_wpre(filepaths: Any, params_dict: Dict, t, x, wpre, u, Ca) -> None:
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
