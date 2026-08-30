"""LFP, rasters, per-pop state traces. standard_plot is the main one.

standard_plot draws the 80/20 weighted-mean LFP, both spike rasters, and
optionally schedule overlays for x0/CE/g_intra/g_inter when those vals are
provided.
"""
import matplotlib.pyplot as plt
import os
import numpy as np
from brian2 import *
from typing import Dict, Any
from config_structs import ZoomConfig

import plotting.style  # noqa: F401
from plotting.plotting_utils import *

def _apply_zoom(axes: list, zoom: ZoomConfig) -> None:
    """Set the x-axis limits on all given axes.

    Args:
        axes: list of matplotlib Axes.
        zoom: ZoomConfig with start and end in seconds.
    
    Returns: void    
    """
    for ax in axes:
        ax.set_xlim(zoom.start, zoom.end)

def find_clim(spike_matrix: np.ndarray) -> float:
    """Max spike count across the matrix, used as the colorbar upper bound.

    Args:
        spike_matrix (np.ndarray): 2D spike-count matrix.

    Returns:
        float: The max value across the matrix.
    """
    return np.max(spike_matrix)

def _vals_to_time_points(vals, t):
    """Map a timed-array value schedule to step-function coordinates over t.

    For each value in vals, generates a [bin_start, bin_end] time segment paired
    with [v, v], then concatenates all segments into a full sparse step representation.

    Args:
        vals: List/array of values.
        t: Time vector (only t[0] and t[-1] are used to set the span).

    Returns:
        Tuple[np.ndarray, np.ndarray]: (time_points, value_points) suitable
        for ax.plot.
    """
    n = len(vals)
    bin_size = (t[-1] - t[0]) / n
    time_edges = t[0] + np.arange(n + 1) * bin_size
    time_segments = [np.array([time_edges[i], time_edges[i+1]]) for i in range(n)]
    val_segments  = [np.array([v, v]) for v in vals]
    return np.concatenate(time_segments), np.concatenate(val_segments)

def standard_plot(filepaths: Any, params_dict: Dict, t: np.ndarray,
        x1: np.ndarray, x2: np.ndarray, spike_matrix_1: np.ndarray, spike_matrix_2: np.ndarray,
        num_cells: int, sim_duration: float, zoom: ZoomConfig = None,
        x0_t=None, ce_t=None,
        g_intra_vals=None, g_inter_vals=None, show: bool = True, png = False) -> None:
    """Weighted-mean LFP + both spike rasters. Adds schedule panels if vals are given.

    Args:
        filepaths (Any): FilePaths with figures_dir.
        params_dict (Dict): Needs SIM_DURATION and TRANSIENT.
        t (np.ndarray): Time vector.
        x1 (np.ndarray): (num_cells, time) HR x.
        x2 (np.ndarray): (num_cells, time) ML x.
        spike_matrix_1 (np.ndarray): HR spike count matrix.
        spike_matrix_2 (np.ndarray): ML spike count matrix.
        num_cells (int): Number of neurons per population.
        sim_duration (float): Simulation duration in seconds.
        zoom: optional ZoomConfig; if provided, restricts x-axis to [start, end].
        x0_t: Optional x0 schedule; adds an extra panel if paired with coupling_vals.
        ce_t: Optional CE schedule; pairs with x_naught_vals.
        g_intra_vals: Optional intra-pop conductance schedule.
        g_inter_vals: Optional inter-pop conductance schedule; pairs with g_intra_vals.
        show (bool): If True, call plt.show() before closing the figure.
        png (bool): If True, saves as a rasterized png. Otherwise, saves as svg.
    
    Returns: void
    """

    has_x0_ce = x0_t is not None and ce_t is not None
    has_g = g_intra_vals is not None and g_inter_vals is not None
    n_axes = 3 + has_x0_ce + has_g
    height_ratios = [3, 3, 3] + ([1] if has_x0_ce else []) + ([1] if has_g else [])

    fig, axes = plt.subplots(n_axes, 1, figsize=(30, 15), sharex=True,
                             layout='constrained',
                             gridspec_kw={'height_ratios': height_ratios})
    ax1, ax2, ax3 = axes[0], axes[1], axes[2]
    ax4 = axes[3] if has_x0_ce else None
    ax5 = axes[3 + has_x0_ce] if has_g else None

    fig.suptitle('Weighted LFP & Spike Rasters')
    fig.set_constrained_layout_pads(w_pad=0.1, h_pad=0.1,
                                     wspace=0.02, hspace=0.02)

    x1_mean = bin_data(np.mean(x1, axis=0), 100)
    x2_mean = bin_data(np.mean(x2, axis=0), 100)
    x_mean = (0.8 * x1_mean) + (0.2 * x2_mean)
    ax1.plot(t, x_mean)
    ax1.set_ylabel("LFP (a.u.)", fontsize=10)
    ax1.set_title("LFP signal (80/20 weight excitatory/inhibitory)", fontsize=10)
    ax1.tick_params(axis='y', labelsize=10)

    # raster 1
    HR_CLIM = max(5, find_clim(spike_matrix_1))
    print(spike_matrix_1.shape)
    raster1 = ax2.imshow(spike_matrix_1, interpolation='none', aspect='auto',
                   origin='lower', extent=[params_dict['TRANSIENT'], sim_duration+params_dict['TRANSIENT'], 0, num_cells], clim=(0, HR_CLIM))

    ax2.set_ylabel('Neuron index', fontsize=10)
    ax2.set_title('Excitatory Population Spike Raster (Spike Count)', fontsize=10)
    ax2.tick_params(axis='y', labelsize=10)
    
    # config colorbar
    cbar = fig.colorbar(raster1, ax=ax2, location='right', aspect=25, pad=0.001)
    cbar.ax.tick_params(labelsize=10)
    cbar.ax.set_ylim(0, HR_CLIM)
    cbar.minorticks_on()

    # raster 2
    ML_CLIM = find_clim(spike_matrix_2)
    raster2 = ax3.imshow(spike_matrix_2, interpolation='none', aspect='auto',
                   origin='lower', extent=[params_dict['TRANSIENT'], sim_duration+params_dict['TRANSIENT'], 0, num_cells], clim=(0, ML_CLIM))

    ax3.set_ylabel('Neuron index', fontsize=10)
    ax3.set_title('Inhibitory Population Spike Raster (Spike Count)', fontsize=10)
    ax3.tick_params(axis='y', labelsize=10)

    # config colorbar
    cbar = fig.colorbar(raster2, ax=ax3, location='right', aspect=25, pad=0.001)
    cbar.ax.tick_params(labelsize=10)
    cbar.minorticks_on()

    if has_x0_ce:
        ax4.plot(t, x0_t, label='x0', color='blue')
        ax4b = ax4.twinx()
        ax4b.plot(t, ce_t, label='Ce', color='orange')
        ax4.set_ylabel("Epop excite", fontsize=10, c='blue')
        ax4b.set_ylabel("Coupling (J)", fontsize=10, c='orange')
        ax4.tick_params(axis='y', labelcolor='blue', labelsize=10)
        ax4b.tick_params(axis='y', labelcolor='orange', labelsize=10)
        ax4.set_title("J and Epop excite over time", fontsize=10)

    if has_g:
        gintra_t, gintra_v = _vals_to_time_points(g_intra_vals, t)
        ginter_t, ginter_v = _vals_to_time_points(g_inter_vals, t)
        ax5.plot(gintra_t, gintra_v, label='g_intra', color='blue')
        ax5.plot(ginter_t, ginter_v, label='g_inter', color='orange')
        ax5.set_title("g variables over time", fontsize=10)
        ax5.set_ylabel("Conductance (uS)", fontsize=10)
        ax5.tick_params(axis='y', labelsize=10)
        ax5.yaxis.get_offset_text().set_fontsize(10)
        ax5.legend()

    axes[len(axes) - 1].set_xlabel("Time (s)", fontsize=10)
    axes[len(axes) - 1].tick_params(axis='x', labelsize=10)

    plt.xlim(params_dict['TRANSIENT'], params_dict['SIM_DURATION']/second + params_dict['TRANSIENT'])
    if zoom:
        _apply_zoom(axes, zoom)
    file_extension = "png" if png else "svg"
    fig.get_layout_engine().set(w_pad=0.2, h_pad=0.2, hspace=0.2, wspace=0.2)
    out_path = os.path.join(filepaths.figures_dir, "standard_plot." + file_extension)
    plt.savefig(out_path, format=file_extension, dpi=300)

    if show:
        plt.show()
    plt.close()

def raster_plot(filepaths: Any, population: int, t: np.ndarray,
        x: np.ndarray, spike_matrix: np.ndarray, num_cells: int, sim_duration: float,
        zoom: ZoomConfig = None) -> None:
    """Mean x trace above the spike raster for one population (1=EPOP, else IPOP).
        
    Args:
        filepaths: FilePaths with figures_dir.
        population: 1 for excitatory, anything else for inhibitory.
        t: time vector.
        x: (num_cells, time) state variable.
        spike_matrix: spike count matrix.
        num_cells: number of neurons.
        sim_duration: simulation duration in seconds.
        zoom: optional ZoomConfig; if provided, restricts x-axis to [start, end].
        
    Returns: void
    """
    population_name = "Excitatory Population" if population == 1 else "Inhibitory Population"
    clim_max = find_clim(spike_matrix)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(30, 9), sharex=True, constrained_layout=True)
    fig.suptitle(f'{population_name} - All Neurons Averaged')

    x_mean = np.mean(x, axis=0)
    ax1.plot(t, x_mean)
    ax1.set_ylabel("Mean x")

    raster = ax2.imshow(spike_matrix, interpolation='none', aspect='auto',
                   origin='lower', extent=[0, sim_duration, 0, num_cells], clim=(0, clim_max))

    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('Neuron index')
    ax2.set_title(f'{population_name} Spike Raster (Spike Count)')

    cbar = fig.colorbar(raster, ax=ax2, location='right', aspect=25, pad=0.001)
    cbar.minorticks_on()

    if zoom:
        _apply_zoom([ax1, ax2], zoom)

    fig.get_layout_engine().set(w_pad=0.2, h_pad=0.2, hspace=0.2, wspace=0.2)
    plt.savefig(os.path.join(filepaths.figures_dir, f"{population_name}_raster.png"), format='png')
    plt.show()

def plot_epop_multiple(filepaths: Any, t: np.ndarray, x1: np.ndarray,
        zoom: ZoomConfig = None) -> None:
    """x-traces for the first 3 EPOP neurons.
    
    Args:
        filepaths: FilePaths with figures_dir.
        t: time vector.
        x1: (num_cells, time) EPOP x.
        zoom: optional ZoomConfig; if provided, restricts x-axis to [start, end].
    
    Returns: void
    """
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(30, 10), sharex=True)
    fig.suptitle("Excitatory Population - Multiple Neurons")
    ax1.plot(t, x1[0])
    ax1.set_ylabel("Neuron 0 x")
    ax2.plot(t, x1[1])
    ax2.set_ylabel("Neuron 1 x")
    ax3.plot(t, x1[2])
    ax3.set_ylabel("Neuron 2 x")
    ax3.set_xlabel("Time (s)")

    if zoom:
        _apply_zoom([ax1, ax2, ax3], zoom)

    plt.savefig(os.path.join(filepaths.figures_dir, "pop1_multiple_neurons.png"), format="png")
    plt.show()

def plot_epop_single(filepaths: Any, t: np.ndarray, x1: np.ndarray, y1: np.ndarray,
        z1: np.ndarray, I_syn_inter_1: np.ndarray, zoom: ZoomConfig = None) -> None:
    """x, y, z, and I_syn_inter for EPOP neuron 0.

    Args:
        filepaths: FilePaths with figures_dir.
        t: time vector.
        x1: (num_cells, time) HR x.
        y1: (num_cells, time) HR y.
        z1: (num_cells, time) HR z.
        I_syn_inter_1: (num_cells, time) inter-pop synaptic current into EPOP.
        zoom: optional ZoomConfig; if provided, restricts x-axis to [start, end].
    
    Returns: void
    """
    fig, (ax1) = plt.subplots(1, 1, figsize=(20, 10), sharex=True)
    fig.suptitle("Excitatory Population Variables - One Neuron")
    ax1.plot(t, x1[0])
    ax1.set_ylabel("x (Neuron 0)")
    # ax2.plot(t, y1[0])
    # ax2.set_ylabel("Neuron 0 y")
    # ax3.plot(t, z1[0])
    # ax3.set_ylabel("Neuron 0 z")
    # ax4.plot(t, I_syn_inter_1[0])
    # ax4.set_ylabel("Neuron 0 I_{syn, inter}")
    ax1.set_xlabel("Time (s)")

    if zoom:
        _apply_zoom([ax1], zoom)

    plt.savefig(os.path.join(filepaths.figures_dir, "pop1_single_neuron.svg"), format="svg")
    plt.show()

def plot_epop_mean(filepaths: Any, t: np.ndarray, x1: np.ndarray, y1: np.ndarray,
        z1: np.ndarray, zoom: ZoomConfig = None) -> None:
    """Mean x, y, z over all EPOP neurons.

    Args:
        filepaths: FilePaths with figures_dir.
        t: time vector.
        x1: (num_cells, time) HR x.
        y1: (num_cells, time) HR y.
        z1: (num_cells, time) HR z.
        zoom: optional ZoomConfig; if provided, restricts x-axis to [start, end].
    
    Returns: void
    """
    x1_mean = np.mean(x1, axis=0)
    y1_mean = np.mean(y1, axis=0)
    z1_mean = np.mean(z1, axis=0)

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(30, 10), sharex=True)
    fig.suptitle("Excitatory Population Variables - All Neurons Averaged")
    ax1.plot(t, x1_mean)
    ax1.set_ylabel("Mean x1")
    ax2.plot(t, y1_mean)
    ax2.set_ylabel("Mean y1")
    ax3.plot(t, z1_mean)
    ax3.set_ylabel("Mean z1")
    ax3.set_xlabel("Time (s)")

    if zoom:
        _apply_zoom([ax1, ax2, ax3], zoom)

    plt.savefig(os.path.join(filepaths.figures_dir, "pop1_mean_neurons.png"), format="png")
    plt.show()

def plot_ipop_single(filepaths: Any, t: np.ndarray, x2: np.ndarray, n2: np.ndarray,
        zoom: ZoomConfig = None) -> None:
    """x and n for IPOP neuron 0.
    
    Args:
        filepaths: FilePaths with figures_dir.
        t: time vector.
        x2: (num_cells, time) ML x.
        n2: (num_cells, time) ML gating variable n.
        zoom: optional ZoomConfig; if provided, restricts x-axis to [start, end].
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 10), sharex=True)
    fig.suptitle("Inhibitory Population Variables - One Neuron")
    ax1.plot(t, x2[0])
    ax1.set_ylabel("x (Neuron 0)")
    ax2.plot(t, n2[0])
    ax2.set_ylabel("Neuron 0 n")
    ax1.set_xlabel("Time (s)")

    if zoom:
        _apply_zoom([ax1], zoom)

    plt.savefig(os.path.join(filepaths.figures_dir, "pop2_single_neuron.svg"), format="svg")
    plt.show()

def plot_ipop_mean():
    """Not implemented."""
    pass

def plot_both(filepaths: Any, t: np.ndarray, x1: np.ndarray, x2: np.ndarray,
        zoom: ZoomConfig = None) -> None:
    """x of neuron 0 from EPOP and IPOP, stacked.

    Args:
        filepaths: FilePaths with figures_dir.
        t: time vector.
        x1: (num_cells, time) HR x.
        x2: (num_cells, time) ML x.
        zoom: optional ZoomConfig; if provided, restricts x-axis to [start, end].
    
    Returns: void
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.suptitle("One Neuron From Excitatory and Inhibitory Populations")
    ax1.plot(t, x1[0])
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("x1")
    ax2.plot(t, x2[0])
    ax2.set_xlabel("Time (s)")
    ax2.set_ylabel("x2")

    if zoom:
        _apply_zoom([ax1, ax2], zoom)

    plt.savefig(os.path.join(filepaths.figures_dir, "pop1_and_pop2_single.png"), format="png")
    plt.show()


def plot_both_avg(filepaths: Any, t: np.ndarray, x1: np.ndarray, y1: np.ndarray,
        z1: np.ndarray, x2: np.ndarray, n: np.ndarray, zoom: ZoomConfig = None) -> None:
    """Mean x of EPOP and IPOP, stacked. y1, z1, n are unused.

    Args:
        filepaths: FilePaths with figures_dir.
        t: time vector.
        x1: (num_cells, time) HR x.
        y1: unused.
        z1: unused.
        x2: (num_cells, time) ML x.
        n: unused.
        zoom: optional ZoomConfig; if provided, restricts x-axis to [start, end].
    
    Returns: void
    """
    x1_mean = np.mean(x1, axis=0)
    x2_mean = np.mean(x2, axis=0)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    fig.suptitle("Mean Neuron From Excitatory and Inhibitory Populations")
    ax1.plot(t, x1_mean)
    ax1.set_xlabel("Time (s)")
    ax1.set_ylabel("x1")
    ax2.plot(t, x2_mean)
    ax2.set_xlabel("Time (s)")
    ax2.set_ylabel("x2")

    if zoom:
        _apply_zoom([ax1, ax2], zoom)

    plt.savefig(os.path.join(filepaths.figures_dir, "pop1_and_pop2_single.png"), format="png")
    plt.show()

def plot_synapse_currents(filepaths: Any, t: np.ndarray, epop_intra: np.ndarray,
        epop_inter: np.ndarray, ipop_intra: np.ndarray, ipop_inter: np.ndarray,
        bin_size = 1000):
    """Plot mean synaptic currents.
    
    Args:
        filepaths (Any): FilePaths object containing figures directory path.
        epop_intra (np.ndarray): Mean current for EPOP intrapopulation synaptic current.
        epop_inter (np.ndarray): Mean current for EPOP interpopulation synaptic current.
        ipop_intra (np.ndarray): Mean current for IPOP intrapopulation synaptic current.
        ipop_inter (np.ndarray): Mean current for IPOP interpopulation synaptic current.
        bin_size (int): number of dt by which to bin data.
    
    Returns: void
    """
    epop_intra_mean = bin_data(np.mean(epop_intra, axis=0), bin_size)
    epop_inter_mean = bin_data(np.mean(epop_inter, axis=0), bin_size)
    ipop_intra_mean = bin_data(np.mean(ipop_intra, axis=0), bin_size)
    ipop_inter_mean = bin_data(np.mean(ipop_inter, axis=0), bin_size)

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(10, 8), sharex=True)
    fig.suptitle("Mean Synaptic Currents (uS)")
    ax1.plot(t, epop_intra_mean)
    ax1.set_ylabel(u"E \u2192 E", fontsize=10)
    ax1.tick_params(axis='both', which='major', labelsize=10)
    ax1.yaxis.get_offset_text().set_fontsize(10)
    ax2.plot(t, epop_inter_mean)
    ax2.set_ylabel(u"E \u2192 I", fontsize=10)
    ax2.tick_params(axis='both', which='major', labelsize=10)
    ax2.yaxis.get_offset_text().set_fontsize(10)
    ax3.plot(t, ipop_intra_mean)
    ax3.set_ylabel(u"I \u2192 I", fontsize=10)
    ax3.tick_params(axis='both', which='major', labelsize=10)
    ax3.yaxis.get_offset_text().set_fontsize(10)
    ax4.plot(t, ipop_inter_mean)
    ax4.set_ylabel(u"I \u2192 E", fontsize=10)
    ax4.tick_params(axis='both', which='major', labelsize=10)
    ax4.yaxis.get_offset_text().set_fontsize(10)
    ax4.set_xlabel("Time (s)", fontsize=10)

    plt.savefig(os.path.join(filepaths.figures_dir, "mean_syn_currents.png"), format="png")
    plt.show()