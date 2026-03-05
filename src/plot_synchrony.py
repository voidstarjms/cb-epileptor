import numpy as np
import matplotlib.pyplot as plt
import os
import config

FIGURES_DIR = config.FIGURES_DIR

def plot_synchrony(chi_mean, chi_sd,
                   param1_values, param2_values,
                   param1_label, param2_label):

    plot_synchrony_single(chi_mean, param1_values, param2_values,
                          param1_label, param2_label,
                          title=r'Synchrony $\chi$',
                          vmin=0, vmax=1,
                          save_name='synchrony_chi_mean.png')

    sd_max = np.nanmax(chi_sd) if np.nanmax(chi_sd) > 0 else 0.15
    plot_synchrony_single(chi_sd, param1_values, param2_values,
                          param1_label, param2_label,
                          title=r'SD of $\chi$',
                          vmin=0, vmax=sd_max,
                          save_name='synchrony_chi_sd.png')


def plot_synchrony_single(chi_matrix, param1_values, param2_values,
                          param1_label, param2_label,
                          title=r'Synchrony $\chi$',
                          vmin=0, vmax=1,
                          save_name='synchrony_single.png'):
    
    if not os.path.exists(FIGURES_DIR):
        os.makedirs(FIGURES_DIR)

    fig, ax = plt.subplots(1, 1, figsize=(7, 6))

    p1_edges = _make_edges(param1_values)
    p2_edges = _make_edges(param2_values)

    im = ax.pcolormesh(p1_edges, p2_edges, chi_matrix,
                       cmap='YlOrRd_r', vmin=vmin, vmax=vmax, shading='flat')
    ax.set_xlabel(param1_label, fontsize=13)
    ax.set_ylabel(param2_label, fontsize=13)
    ax.set_title(title, fontsize=14, fontweight='bold')
    fig.colorbar(im, ax=ax, pad=0.02)

    fig.tight_layout()
    filepath = os.path.join(FIGURES_DIR, save_name)
    fig.savefig(filepath, dpi=200, bbox_inches='tight')
    plt.show()
    print(f"Saved to {filepath}")


def _make_edges(values):
    """Convert cell-center values to cell-edge values for pcolormesh"""
    values = np.asarray(values, dtype=float)
    if len(values) < 2:
        return np.array([values[0] - 0.5, values[0] + 0.5])
    half = (values[1] - values[0]) / 2.0
    edges = np.concatenate([
        [values[0] - half],
        (values[:-1] + values[1:]) / 2.0,
        [values[-1] + half]
    ])
    return edges