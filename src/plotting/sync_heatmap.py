"""Sweep heatmaps: mean chi and SD chi over a (param1, param2) grid."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import os
from run import FilePaths

import plotting.style  # noqa: F401

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
    ax.set_xlabel(param1_label)
    ax.set_ylabel(param2_label)
    ax.set_title(title, fontweight='bold')
    fig.colorbar(im, ax=ax, pad=0.02)

    fig.tight_layout()
    filepath = os.path.join(filepaths.figures_dir, save_name)
    fig.savefig(filepath, dpi=200, bbox_inches='tight')
    plt.show()
    print(f"Saved to {filepath}")

def plot_synchrony_multifaceted(filepaths, chi_mean, chi_sd,
                                param1_values, param2_values, param3_values,
                                param4_values, param1_label, param2_label,
                                param3_label, param4_label, run_num=1, plot_mode='r'):
    """Write the multifaceted mean-chi and SD-chi heatmaps, both prefixed with run_num.
        INPUT:
            filepaths: FilePaths with figures_dir.
            chi_mean: 4D array of mean chi values, shape (len(param2), len(param1),
                len(param4), len(param3)).
            chi_sd: 4D array of chi SDs, same shape as chi_mean.
            param1_values: values along the outer x axis.
            param2_values: values along the outer y axis.
            param3_values: values along the inner x axis.
            param4_values: values along the inner y axis.
            param1_label: outer x-axis label (LaTeX allowed).
            param2_label: outer y-axis label (LaTeX allowed).
            param3_label: inner x-axis label (LaTeX allowed).
            param4_label: inner y-axis label (LaTeX allowed).
            run_num: integer prefix for the output filenames.
    """
    save_name_prefix = f'{run_num}_synchrony_r' if plot_mode == 'r' else f'{run_num}_synchrony_chi'
    plot_sync_metric_name = 'r' if plot_mode == 'r' else r'$\chi$'
    plot_synchrony_multifaceted_single(filepaths, chi_mean, param1_values, param2_values,
                          param3_values, param4_values, param1_label, param2_label,
                          param3_label, param4_label, title='Synchrony '+plot_sync_metric_name,
                          vmin=0, vmax=1,
                          save_name=save_name_prefix+'_mean_full.png')

    sd_max = np.nanmax(chi_sd) if np.nanmax(chi_sd) > 0 else 0.15
    plot_synchrony_multifaceted_single(filepaths, chi_sd, param1_values, param2_values,
                          param3_values, param4_values, param1_label, param2_label,
                          param3_label, param4_label, title='SD of '+plot_sync_metric_name,
                          vmin=0, vmax=sd_max,
                          save_name=save_name_prefix+'_sd_full.png')


def plot_synchrony_multifaceted_single(filepaths, chi_matrix, param1_values, param2_values,
                                param3_values, param4_values, param1_label,
                                param2_label, param3_label, param4_label,
                                title, vmin=0, vmax=1,
                                save_name='synchrony_all_params.png'):

    n_outer_rows = len(param1_values)
    n_outer_cols = len(param2_values)

    OUTER_LEFT = 0.18
    OUTER_BOTTOM = 0.05
    OUTER_RIGHT = 0.1
    OUTER_TOP = 0.12

    HSPACE = 0.15
    WSPACE = 0.15

    OUTER_AXIS_OFFSET = 0.06

    fig = plt.figure(figsize=(16, 13))

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
    print(param1_values)
    for i, p1 in enumerate(param1_values):
        for j, p2 in enumerate(param2_values):
            ax = axes[i, j]
            data = chi_matrix[i, j, :, :]

            im = ax.imshow(
                data.T,
                origin="lower",
                aspect="auto",
                cmap='nipy_spectral',
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
                ax.set_xlabel(param3_label, fontsize=10, labelpad=4)
            if j == 0:
                ax.set_ylabel(param4_label, fontsize=10, labelpad=4)

            for spine in ax.spines.values():
                spine.set_linewidth(0.5)

    gs_left = gs.left
    gs_right = gs.right
    gs_top = gs.top
    gs_bottom = gs.bottom

    row_height = (gs_top - gs_bottom) / n_outer_rows

    for i, p1 in enumerate(param1_values):
        row_center = gs_bottom + (i + 0.5) * row_height

        fig.add_artist(plt.Line2D(
            [gs_left - OUTER_AXIS_OFFSET - 0.01, gs_left - OUTER_AXIS_OFFSET],
            [row_center, row_center],
            transform=fig.transFigure,
            color="black", linewidth=0.8, clip_on=False,
        ))

        fig.text(
            gs_left - OUTER_AXIS_OFFSET - 0.015, row_center,
            str(p1),
            ha="right", va="center", fontsize=9,
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
            str(p2),
            ha="center", va="bottom", fontsize=9,
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
        gs_left - 0.15, (gs_top + gs_bottom) / 2,
        param1_label,
        ha="center", va="top", fontsize=11, fontweight="bold",
        transform=fig.transFigure, rotation=90,
    )
    
    # Param 2 label
    fig.text(
        (gs_left + gs_right) / 2, gs_top + 0.05,
        param2_label,
        ha="center", va="center", fontsize=11, fontweight="bold",
        transform=fig.transFigure,
    )

    cbar_ax = fig.add_axes([
        1 - OUTER_RIGHT + 0.01,
        gs_bottom,
        0.025,
        gs_top - gs_bottom,
    ])
    cb = fig.colorbar(ScalarMappable(cmap='nipy_spectral'), cax=cbar_ax)
    cb.set_label("Synchrony", fontsize=11, fontweight="bold")
    cb.ax.tick_params(labelsize=8)

    fig.suptitle(title, fontsize=13, fontweight="bold", y=0.98)

    fname = os.path.join(filepaths.figures_dir, save_name)
    plt.savefig(fname, dpi=200, bbox_inches="tight")
    plt.show()


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
