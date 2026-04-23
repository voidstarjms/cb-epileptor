import matplotlib.pyplot as plt
import os
import numpy as np
from brian2 import *
from typing import Dict, Any

import plotting.style  # noqa: F401


def _apply_zoom(axes: list, params_dict: Dict) -> None:
    XMIN, XMAX = params_dict['SIM_DURATION']/second - 15, params_dict['SIM_DURATION']/second - 10
    for ax in axes:
        ax.set_xlim(XMIN, XMAX)

def find_clim(spike_matrix: np.ndarray) -> float:
    return np.max(spike_matrix)

def _vals_to_time_points(vals, t):
    """Map a timed array value schedule to step-function coordinates over t.
    For each value, generates a [bin_start, bin_end] time segment paired with
    [v, v], then concatenates all segments into a full sparse step representation."""
    n = len(vals)
    bin_size = (t[-1] - t[0]) / n
    time_edges = t[0] + np.arange(n + 1) * bin_size
    time_segments = [np.array([time_edges[i], time_edges[i+1]]) for i in range(n)]
    val_segments  = [np.array([v, v]) for v in vals]
    return np.concatenate(time_segments), np.concatenate(val_segments)

def standard_plot(filepaths: Any, params_dict: Dict, t: np.ndarray,
        x1: np.ndarray, x2: np.ndarray, spike_matrix_1: np.ndarray, spike_matrix_2: np.ndarray,
        num_cells: int, sim_duration: float, zoom: bool = False,
        x_naught_vals=None, coupling_vals=None,
        g_intra_vals=None, g_inter_vals=None) -> None:

    has_x0_ce = x_naught_vals is not None and coupling_vals is not None
    has_g = g_intra_vals is not None and g_inter_vals is not None
    n_axes = 3 + has_x0_ce + has_g
    height_ratios = [3, 3, 3] + ([1] if has_x0_ce else []) + ([1] if has_g else [])

    fig, axes = plt.subplots(n_axes, 1, figsize=(30, 15), sharex=True,
                             layout='constrained',
                             gridspec_kw={'height_ratios': height_ratios})
    ax1, ax2, ax3 = axes[0], axes[1], axes[2]
    ax4 = axes[3] if has_x0_ce else None
    ax5 = axes[3 + has_x0_ce] if has_g else None

    fig.suptitle('Weighted LFP + Both Rasters')
    fig.set_constrained_layout_pads(w_pad=0.1, h_pad=0.1,
                                     wspace=0.02, hspace=0.02)

    x1_mean = np.mean(x1, axis=0)
    x2_mean = np.mean(x2, axis=0)
    x_mean = (0.8 * x1_mean) + (0.2 * x2_mean)
    ax1.plot(t, x_mean)
    ax1.set_ylabel("Mean x weighted 80/20")
    ax1.set_title("LFP signal (80/20 weight excitatory/inhibitory)")

    # raster 1
    HR_CLIM = find_clim(spike_matrix_1)
    raster1 = ax2.imshow(spike_matrix_1, interpolation='none', aspect='auto',
                   origin='lower', extent=[params_dict['TRANSIENT'], sim_duration+params_dict['TRANSIENT'], 0, num_cells], clim=(0, HR_CLIM))

    ax2.set_ylabel('Neuron index')
    ax2.set_title('Excitatory Population Spike Raster (Spike Count)')

    # config colorbar
    cbar = fig.colorbar(raster1, ax=ax2, location='right', aspect=25, pad=0.001)
    cbar.minorticks_on()

    # raster 2
    ML_CLIM = find_clim(spike_matrix_2)
    raster2 = ax3.imshow(spike_matrix_2, interpolation='none', aspect='auto',
                   origin='lower', extent=[params_dict['TRANSIENT'], sim_duration+params_dict['TRANSIENT'], 0, num_cells], clim=(0, ML_CLIM))

    # ax3.set_xlabel('Time (s)', fontsize=12)
    ax3.set_ylabel('Neuron index')
    ax3.set_title('Inhibitory Population Spike Raster (Spike Count)')

    # config colorbar
    cbar = fig.colorbar(raster2, ax=ax3, location='right', aspect=25, pad=0.001)
    cbar.minorticks_on()

    # plot timed arrays if they exist
    if has_x0_ce:
        x0_t, x0_v = _vals_to_time_points(x_naught_vals, t)
        ce_t, ce_v = _vals_to_time_points(coupling_vals, t)
        ax4.plot(x0_t, x0_v, label='x0', color='blue')
        ax4.plot(ce_t, ce_v, label='Ce', color='orange')
        ax4.set_ylabel("x0")
        ax4.set_title("Ce and x0 over time")
        ax4.legend()

    if has_g:
        gintra_t, gintra_v = _vals_to_time_points(g_intra_vals, t)
        ginter_t, ginter_v = _vals_to_time_points(g_inter_vals, t)
        ax5.plot(gintra_t, gintra_v, label='g_intra', color='blue')
        ax5.plot(ginter_t, ginter_v, label='g_inter', color='orange')
        ax5.set_title("g variables over time")
        ax5.set_xlabel("Time (s)")
        ax5.set_ylabel("Conductance (uS)")
        ax5.legend()

    plt.xlim(params_dict['TRANSIENT'], params_dict['SIM_DURATION']/second + params_dict['TRANSIENT'])
    # optionally zoom
    if zoom:
        _apply_zoom(axes, params_dict)
    # save plot
    fig.get_layout_engine().set(w_pad=0.2, h_pad=0.2, hspace=0.2, wspace=0.2)
    out_path = os.path.join(filepaths.figures_dir, "standard_plot.png")
    plt.savefig(out_path, format='png')

    plt.show()
    plt.close()

def raster_plot(filepaths: Any, params_dict: Dict, population: int, t: np.ndarray,
        x: np.ndarray, spike_matrix: np.ndarray, num_cells: int, sim_duration: float, zoom: bool = False) -> None:
    population_name = ""
    if population == 1:
        population_name = "Excitatory Population"
    else:
        population_name = "Inhibitory Population"
    clim_max = find_clim(spike_matrix)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(30, 9), sharex=True, constrained_layout=True)
    fig.suptitle(f'{population_name} - All Neurons Averaged')

    # Plot the averaged data instead of just neuron 0
    x_mean = np.mean(x, axis=0)
    ax1.plot(t, x_mean)
    ax1.set_ylabel("Mean x")

    # configure main raster plot
    raster = ax2.imshow(spike_matrix, interpolation='none', aspect='auto',
                   origin='lower', extent=[0, sim_duration, 0, num_cells], clim=(0, clim_max))

    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('Neuron index')
    ax2.set_title(f'{population_name} Spike Raster (Spike Count)')

    # config colorbar
    cbar = fig.colorbar(raster, ax=ax2, location='right', aspect=25, pad=0.001)
    cbar.minorticks_on()

    # optionally zoom
    if zoom:
        _apply_zoom([ax1, ax2], params_dict)

    # save plot
    fig.get_layout_engine().set(w_pad=0.2, h_pad=0.2, hspace=0.2, wspace=0.2)
    plt.savefig(os.path.join(filepaths.figures_dir, f"{population_name}_raster.png"), format='png')
    plt.show()

def plot_hr_multiple(filepaths: Any, t: np.ndarray, x1: np.ndarray) -> None:
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(30, 10), sharex=True)
    fig.suptitle("Excitatory Population - Multiple Neurons")
    ax1.plot(t, x1[0])
    ax1.set_ylabel("Neuron 0 x")
    ax2.plot(t, x1[1])
    ax2.set_ylabel("Neuron 1 x")
    ax3.plot(t, x1[2])
    ax3.set_ylabel("Neuron 2 x")
    ax1.legend()

    plt.savefig(os.path.join(filepaths.figures_dir, "pop1_multiple_neurons.png"), format="png")
    plt.show()

def plot_hr_single(filepaths: Any, t: np.ndarray, x1: np.ndarray, y1: np.ndarray,
        z1: np.ndarray, I_syn_inter_1: np.ndarray) -> None:
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(30, 10), sharex=True)
    fig.suptitle("Excitatory Population Variables - One Neuron")
    ax1.plot(t, x1[0])
    ax1.set_ylabel("Neuron 0 x")
    ax2.plot(t, y1[0])
    ax2.set_ylabel("Neuron 0 y")
    ax3.plot(t, z1[0])
    ax3.set_ylabel("Neuron 0 z")
    ax4.plot(t, I_syn_inter_1[0])
    ax4.set_ylabel("Neuron 0 I_{syn, inter}")
    ax4.set_xlabel("Time (s)")
    plt.savefig(os.path.join(filepaths.figures_dir, "pop1_single_neuron.png"), format="png")
    plt.show()

def plot_hr_mean(filepaths: Any, t: np.ndarray, x1: np.ndarray, y1: np.ndarray, z1: np.ndarray) -> None:
    # Calculate the mean across all neurons (axis=0)
    x1_mean = np.mean(x1, axis=0)
    y1_mean = np.mean(y1, axis=0)
    z1_mean = np.mean(z1, axis=0)

    # All pop 1 variables (Figure 2 - Now Averaged)
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(30, 10), sharex=True)
    fig.suptitle("Excitatory Population Variables - All Neurons Averaged")

    # Plot the averaged data instead of just neuron 0
    ax1.plot(t, x1_mean)
    ax1.set_ylabel("Mean x1")
    ax2.plot(t, y1_mean)
    ax2.set_ylabel("Mean y1")
    ax3.plot(t, z1_mean)
    ax3.set_ylabel("Mean z1")
    ax3.set_xlabel("Time (s)")

    plt.savefig(os.path.join(filepaths.figures_dir, "pop1_mean_neurons.png"), format="png")
    plt.show()

def plot_ml_single(filepaths: Any, t: np.ndarray, x2: np.ndarray, n2: np.ndarray) -> None:
    # All pop2 variables
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(30, 10), sharex=True)
    fig.suptitle("Inhibitory Population Variables - One Neuron")
    ax1.plot(t, x2[0])
    ax1.set_ylabel("Neuron 0 x")
    ax2.plot(t, n2[0])
    ax2.set_ylabel("Neuron 0 n")
    ax2.set_xlabel("Time (s)")
    plt.savefig(os.path.join(filepaths.figures_dir, "pop2_single_neuron.png"), format="png")
    plt.show()

def plot_ml_mean():
    pass

def plot_both(filepaths: Any, t: np.ndarray, x1: np.ndarray, x2: np.ndarray) -> None:
    # One neuron from both pops
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.suptitle("One Neuron From Excitatory and Inhibitory Populations")
    ax1.plot(t, x1[0])
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("x1")
    ax2.plot(t, x2[0])
    ax2.set_xlabel("Time (s)")
    ax2.set_ylabel("x2")

    plt.savefig(os.path.join(filepaths.figures_dir, "pop1_and_pop2_single.png"), format="png")
    plt.show()


def plot_both_avg(filepaths: Any, t: np.ndarray, x1: np.ndarray, y1: np.ndarray,
        z1: np.ndarray, x2: np.ndarray, n: np.ndarray) -> None:
    # neurons averaged from both pops
    x1_mean = np.mean(x1, axis=0)
    y1_mean = np.mean(y1, axis=0)
    z1_mean = np.mean(z1, axis=0)
    x2_mean = np.mean(x2, axis=0)
    n_mean = np.mean(n, axis=0)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.suptitle("Mean Neuron From Excitatory and Inhibitory Populations")
    ax1.plot(t, x1_mean)
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("x1")
    ax2.plot(t, x2_mean)
    ax2.set_xlabel("Time (s)")
    ax2.set_ylabel("x2")

    plt.savefig(os.path.join(filepaths.figures_dir, "pop1_and_pop2_single.png"), format="png")
    plt.show()
