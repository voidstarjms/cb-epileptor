"""Sweep heatmaps for sweep/plot_sweep.py.

plot_synchrony_single draws one 2D grid; plot_synchrony_multifaceted_single
draws a grid of 2D panels for 4 swept parameters. Callers pass the title,
colour range and colorbar label, so nothing here assumes a particular metric.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import os

import plotting.style  # noqa: F401

def plot_synchrony_single(filepaths, chi_matrix, param1_values, param2_values,
                          param1_label, param2_label,
                          title=r'Synchrony $\chi$',
                          vmin=0, vmax=1,
                          save_name='synchrony_single.png', cbar_label=''):
    """pcolormesh of chi_matrix with param1 on x, param2 on y.

    Args:
        filepaths: FilePaths with figures_dir.
        chi_matrix (np.ndarray): 2D array, shape (len(param2), len(param1)).
        param1_values: x-axis values.
        param2_values: y-axis values.
        param1_label (str): x-axis label.
        param2_label (str): y-axis label.
        title (str): Plot title.
        vmin (float): Colormap lower bound.
        vmax (float): Colormap upper bound.
        save_name (str): Output filename inside filepaths.figures_dir.
        cbar_label (str): Colorbar label.
    """
    if not os.path.exists(filepaths.figures_dir):
        os.makedirs(filepaths.figures_dir)

    fig, ax = plt.subplots(1, 1, figsize=(7, 6))

    p1_edges = _make_edges(param1_values)
    p2_edges = _make_edges(param2_values)

    im = ax.pcolormesh(p1_edges, p2_edges, chi_matrix,
                       cmap='viridis', vmin=vmin, vmax=vmax, shading='flat')
    # Ticks on the swept values (cell centres), not matplotlib's round numbers.
    ax.set_xticks(param1_values, [f'{v:.4g}' for v in param1_values])
    ax.set_yticks(param2_values, [f'{v:.4g}' for v in param2_values])
    ax.set_xlabel(param1_label)
    ax.set_ylabel(param2_label)
    ax.set_title(title, fontweight='bold')
    fig.colorbar(im, ax=ax, pad=0.02, label=cbar_label)

    fig.tight_layout()
    filepath = os.path.join(filepaths.figures_dir, save_name)
    fig.savefig(filepath, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved to {filepath}")

def plot_synchrony_multifaceted_single(filepaths, chi_matrix, param1_values, param2_values,
                                param3_values, param4_values, param1_label,
                                param2_label, param3_label, param4_label,
                                title, vmin=0, vmax=1,
                                save_name='synchrony_all_params.png', cbar_label=''):
    """Grid of heatmap panels: param1 picks the panel row, param2 the panel column,
    and each panel shows param3 (x) against param4 (y).

    Args:
        filepaths: Object with figures_dir.
        chi_matrix (np.ndarray): 4D array indexed [param1, param2, param3, param4].
        param1_values .. param4_values: Values along each of those axes.
        param1_label .. param4_label (str): Axis labels (LaTeX allowed).
        title (str): Figure title.
        vmin, vmax (float): Colour range shared by every panel.
        save_name (str): Output filename inside filepaths.figures_dir.
        cbar_label (str): Colorbar label.
    """
    os.makedirs(filepaths.figures_dir, exist_ok=True)

    n_outer_rows = len(param1_values)
    n_outer_cols = len(param2_values)

    OUTER_LEFT = 0.13
    OUTER_BOTTOM = 0.07
    OUTER_RIGHT = 0.1
    OUTER_TOP = 0.13

    HSPACE = 0.15
    WSPACE = 0.15

    OUTER_AXIS_OFFSET = 0.045

    fig = plt.figure(figsize=(15, 15))

    gs = gridspec.GridSpec(
        n_outer_rows, n_outer_cols,
        figure=fig,
        left=OUTER_LEFT, right=1 - OUTER_RIGHT,
        bottom=OUTER_BOTTOM, top=1 - OUTER_TOP,
        hspace=HSPACE, wspace=WSPACE,
    )

    norm = Normalize(vmin=vmin, vmax=vmax)

    axes = np.empty((n_outer_rows, n_outer_cols), dtype=object)
    for i in range(n_outer_rows):
        for j in range(n_outer_cols):
            axes[i, j] = fig.add_subplot(gs[i, j])
    for i, p1 in enumerate(param1_values):
        for j, p2 in enumerate(param2_values):
            ax = axes[i, j]
            data = chi_matrix[i, j, :, :]
            # print(data)

            im = ax.imshow(
                data.T,
                origin="lower",
                aspect="auto",
                cmap='viridis',
                norm=norm,
            )

            ax.set_xticks(range(len(param3_values)))
            ax.set_xticklabels(
                [str(v) for v in param3_values],
                fontsize=7, rotation=45, ha="right",
            )

            ax.set_yticks(range(len(param4_values)))
            ax.set_yticklabels([str(v) for v in param4_values], fontsize=7)

            if i == n_outer_rows - 1:
                ax.set_xlabel(param3_label, fontsize=15, labelpad=4)
            if j == 0:
                ax.set_ylabel(param4_label, fontsize=15, labelpad=4)

            for spine in ax.spines.values():
                spine.set_linewidth(0.5)

    gs_left = gs.left
    gs_right = gs.right
    gs_top = gs.top
    gs_bottom = gs.bottom

    row_height = (gs_top - gs_bottom) / n_outer_rows

    for i, p1 in enumerate(param1_values):
        row_center = gs_top - (i + 0.5) * row_height

        fig.add_artist(plt.Line2D(
            [gs_left - OUTER_AXIS_OFFSET - 0.01, gs_left - OUTER_AXIS_OFFSET],
            [row_center, row_center],
            transform=fig.transFigure,
            color="black", linewidth=0.8, clip_on=False,
        ))

        fig.text(
            gs_left - OUTER_AXIS_OFFSET - 0.015, row_center,
            '%.3f' % p1,
            ha="right", va="center", fontsize=11,
            transform=fig.transFigure,
        )

    fig.add_artist(plt.Line2D(
        [gs_left - OUTER_AXIS_OFFSET, gs_left - OUTER_AXIS_OFFSET],
        [gs_bottom, gs_top],
        transform=fig.transFigure,
        color="black", linewidth=0.8, clip_on=False,
    ))

    col_width = (gs_right - gs_left) / n_outer_cols

    for j, p2 in enumerate(param2_values):
        col_center = gs_left + (j + 0.5) * col_width

        fig.add_artist(plt.Line2D(
            [col_center, col_center],
            [gs_top + 0.005, gs_top + 0.015],
            transform=fig.transFigure,
            color="black", linewidth=0.8, clip_on=False,
        ))

        fig.text(
            col_center, gs_top + 0.025,
            '%.3f' % p2,
            ha="center", va="bottom", fontsize=11,
            transform=fig.transFigure,
        )

    fig.add_artist(plt.Line2D(
        [gs_left, gs_right],
        [gs_top + 0.015, gs_top + 0.015],
        transform=fig.transFigure,
        color="black", linewidth=0.8, clip_on=False,
    ))
    
    # Param 1 label
    fig.text(
        gs_left - 0.12, (gs_top + gs_bottom) / 2,
        param1_label,
        ha="left", va="center", fontsize=11, fontweight="bold",
        transform=fig.transFigure, rotation=90,
    )
    
    # Param 2 label
    fig.text(
        (gs_left + gs_right) / 2, gs_top + 0.07,
        param2_label,
        ha="center", va="top", fontsize=11, fontweight="bold",
        transform=fig.transFigure,
    )

    cbar_ax = fig.add_axes([
        1 - OUTER_RIGHT + 0.01,
        gs_bottom,
        0.025,
        gs_top - gs_bottom,
    ])
    cbar_ax.set_ylim(vmin, vmax)
    sm = ScalarMappable(cmap='viridis')
    sm.set_clim(vmin, vmax)
    cb = fig.colorbar(sm, cax=cbar_ax)
    cb.set_label(cbar_label, fontsize=11, fontweight="bold")
    cb.ax.tick_params(labelsize=8)

    fig.suptitle(title, fontsize=24, fontweight="bold", y=0.98)

    fname = os.path.join(filepaths.figures_dir, save_name)
    fig.savefig(fname, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved to {fname}")


def _make_edges(values):
    """Convert cell-center values to cell-edge values for pcolormesh.

    Args:
        values (np.ndarray): 1D array of cell centers.

    Returns:
        np.ndarray: 1D array of cell edges, length len(values)+1.
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
