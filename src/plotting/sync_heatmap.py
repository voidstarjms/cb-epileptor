"""Sweep heatmaps: mean chi and SD chi over a (param1, param2) grid."""
import numpy as np
import matplotlib.pyplot as plt
import os

def plot_synchrony(filepaths, chi_mean, chi_sd,
                   param1_values, param2_values,
                   param1_label, param2_label,
                   run_num=1):
    """Write the mean-chi and SD-chi heatmaps, both prefixed with run_num.
        INPUT:
            filepaths: FilePaths with figures_dir.
            chi_mean: 2D array of mean chi values, shape (len(param2), len(param1)).
            chi_sd: 2D array of chi SDs, same shape as chi_mean.
            param1_values: values along the x axis.
            param2_values: values along the y axis.
            param1_label: x-axis label (LaTeX allowed).
            param2_label: y-axis label (LaTeX allowed).
            run_num: integer prefix for the output filenames.
    """
    plot_synchrony_single(filepaths, chi_mean, param1_values, param2_values,
                          param1_label, param2_label,
                          title=r'Synchrony $\chi$',
                          vmin=0, vmax=1,
                          save_name=f'{run_num}_synchrony_chi_mean.png')

    sd_max = np.nanmax(chi_sd) if np.nanmax(chi_sd) > 0 else 0.15
    plot_synchrony_single(filepaths, chi_sd, param1_values, param2_values,
                          param1_label, param2_label,
                          title=r'SD of $\chi$',
                          vmin=0, vmax=sd_max,
                          save_name=f'{run_num}_synchrony_chi_sd.png')


def plot_synchrony_single(filepaths, chi_matrix, param1_values, param2_values,
                          param1_label, param2_label,
                          title=r'Synchrony $\chi$',
                          vmin=0, vmax=1,
                          save_name='synchrony_single.png'):
    """pcolormesh of chi_matrix with param1 on x, param2 on y.
        INPUT:
            filepaths: FilePaths with figures_dir.
            chi_matrix: 2D array, shape (len(param2), len(param1)).
            param1_values: x-axis values.
            param2_values: y-axis values.
            param1_label: x-axis label.
            param2_label: y-axis label.
            title: plot title.
            vmin: colormap lower bound.
            vmax: colormap upper bound.
            save_name: output filename inside filepaths.figures_dir.
    """
    if not os.path.exists(filepaths.figures_dir):
        os.makedirs(filepaths.figures_dir)

    fig, ax = plt.subplots(1, 1, figsize=(7, 6))

    p1_edges = _make_edges(param1_values)
    p2_edges = _make_edges(param2_values)

    im = ax.pcolormesh(p1_edges, p2_edges, chi_matrix,
                       cmap='YlOrRd', vmin=vmin, vmax=vmax, shading='flat')
    ax.set_xlabel(param1_label, fontsize=13)
    ax.set_ylabel(param2_label, fontsize=13)
    ax.set_title(title, fontsize=14, fontweight='bold')
    fig.colorbar(im, ax=ax, pad=0.02)

    fig.tight_layout()
    filepath = os.path.join(filepaths.figures_dir, save_name)
    fig.savefig(filepath, dpi=200, bbox_inches='tight')
    plt.show()
    print(f"Saved to {filepath}")


def _make_edges(values):
    """Convert cell-center values to cell-edge values for pcolormesh.
        INPUT:
            values: 1D array of cell centers.
        RETURN:
            1D array of cell edges, length len(values)+1.
    """
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
