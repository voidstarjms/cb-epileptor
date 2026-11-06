"""Plot sweep heatmaps from a results table (made by build_table.py).

    python plot_sweep.py <results.csv> [--metric r_exc r_inh ...] [--vmin V] [--vmax V]
                         [--config PATH] [--out-dir DIR]

For each metric, averages the realizations at every grid point and writes a
mean figure and an SD figure. The plot type depends on how many parameters were swept:
    1 -> line plot with SD error bars
    2 -> heatmap: 1st on x, 2nd on y
    4 -> grid of heatmaps: 3rd picks the panel row, 4th the panel column;
         inside each panel the 1st is x and the 2nd is y

Plot order and labels come from sweep_config.yaml (default: next to the CSV).
"""
import argparse
import os
import sys
from typing import List, Sequence, Tuple

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import plotting.sync_heatmap as heatmap
from sweep_config import load_sweep_config, sweep_param_keys

METRIC_LABELS = {
    'r_exc': 'KOP r (E pop)', 'r_inh': 'KOP r (I pop)',
    'chi_exc': r'$\chi$ (E pop)', 'chi_inh': r'$\chi$ (I pop)',
    'rate_exc_hz': 'Firing rate, E pop (Hz)', 'rate_inh_hz': 'Firing rate, I pop (Hz)',
    'median_freq_exc_hz': 'Median frequency, E pop (Hz)',
    'median_freq_inh_hz': 'Median frequency, I pop (Hz)',
    'isi_cv_exc': 'ISI CV (E pop)', 'isi_cv_inh': 'ISI CV (I pop)',
}


def grid_stats(table: pd.DataFrame, keys: List[str], metric: str) -> pd.DataFrame:
    """Mean, SD and count of `metric` at each grid point (over realizations)."""
    return table.groupby(keys)[metric].agg(['mean', 'std', 'count'])


def to_grid(stats: pd.DataFrame, column: str, order: Sequence[str]) -> Tuple[np.ndarray, List[np.ndarray]]:
    """Arrange one column of grid_stats into an array with dimensions in `order`.

    Grid points with no result become NaN. Returns the array and each
    dimension's sorted values, e.g. order=[y, x] gives array[y_index, x_index].
    """
    flat = stats[column].reset_index().set_index(list(order))[column]
    values = [np.sort(flat.index.get_level_values(key).unique().to_numpy()) for key in order]
    if len(order) == 1:
        full_index = pd.Index(values[0], name=order[0])
    else:
        full_index = pd.MultiIndex.from_product(values, names=list(order))
    array = flat.reindex(full_index).to_numpy().reshape([len(v) for v in values])
    return array, values


def colour_range(array: np.ndarray, vmin, vmax) -> Tuple[float, float]:
    """Use the requested limits, else the data's range."""
    lo = float(np.nanmin(array) if vmin is None else vmin)
    hi = float(np.nanmax(array) if vmax is None else vmax)
    if hi <= lo:  # constant grid: give the colour scale some width
        hi = lo + 1e-9
    return lo, hi


def plot_1d(stats, x, metric_label, out_dir, name):
    mean, (x_values,) = to_grid(stats, 'mean', [x['key']])
    sd, _ = to_grid(stats, 'std', [x['key']])
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.errorbar(x_values, mean, yerr=sd, marker='o', capsize=3)
    ax.set_xlabel(x['label'])
    ax.set_ylabel(metric_label)
    ax.set_title(f'{metric_label}: mean ± SD over realizations', fontweight='bold')
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f'{name}.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved to {path}")


def plot_2d(stats, sweep_params, metric_label, out_dir, name, vmin, vmax, with_sd):
    x, y = sweep_params
    figures = [('mean', metric_label, vmin)]
    if with_sd:
        figures.append(('std', f'SD of {metric_label}', 0))
    for stat, title, lo in figures:
        array, (y_values, x_values) = to_grid(stats, stat, [y['key'], x['key']])
        lo, hi = colour_range(array, lo, vmax if stat == 'mean' else None)
        heatmap.plot_heatmap(
            out_dir, array, x_values, y_values,
            x['label'], y['label'], title=title, vmin=lo, vmax=hi,
            save_name=f'{name}_{stat}.png', cbar_label=title)


def plot_4d(stats, sweep_params, metric_label, out_dir, name, vmin, vmax, with_sd):
    x, y, row, col = sweep_params
    figures = [('mean', metric_label, vmin)]
    if with_sd:
        figures.append(('std', f'SD of {metric_label}', 0))
    for stat, title, lo in figures:
        array, (row_values, col_values, x_values, y_values) = to_grid(
            stats, stat, [row['key'], col['key'], x['key'], y['key']])
        lo, hi = colour_range(array, lo, vmax if stat == 'mean' else None)
        heatmap.plot_heatmap_grid(
            out_dir, array,
            row_values, col_values, x_values, y_values,
            row['label'], col['label'], x['label'], y['label'],
            title=title, vmin=lo, vmax=hi,
            save_name=f'{name}_{stat}.png', cbar_label=title)


def plot_metric(table, config, metric, out_dir, vmin, vmax):
    """Write the mean figure, and the SD figure if there are repeat realizations, for one metric."""
    if metric not in table.columns:
        sys.exit(f"--metric {metric}: no such column")
    sweep_params = config['sweep_params']
    keys = sweep_param_keys(config)
    stats = grid_stats(table, keys, metric)
    label = METRIC_LABELS.get(metric, metric)
    name = f"{metric}_{'_'.join(keys)}"
    # SD needs at least 2 realizations: pandas gives NaN for the SD of one value.
    with_sd = stats['count'].max() > 1
    if not with_sd:
        print(f"  {metric}: one realization per grid point, so no SD figure")
    if len(sweep_params) == 1:
        plot_1d(stats, sweep_params[0], label, out_dir, name)
    elif len(sweep_params) == 2:
        plot_2d(stats, sweep_params, label, out_dir, name, vmin, vmax, with_sd)
    elif len(sweep_params) == 4:
        plot_4d(stats, sweep_params, label, out_dir, name, vmin, vmax, with_sd)
    else:
        sys.exit("Plotting 3 swept parameters isn't supported")
    print(f"  {metric}: jobs per grid point min {stats['count'].min()}, max {stats['count'].max()}")


def main(argv=None):
    parser = argparse.ArgumentParser(description='Plot sweep heatmaps from results.csv.')
    parser.add_argument('table', help='results.csv from build_table.py')
    parser.add_argument('--metric', nargs='+', default=['r_exc', 'r_inh'],
                        help='Column(s) to plot (default: r_exc r_inh)')
    parser.add_argument('--vmin', type=float, help='Colour scale minimum for mean plots')
    parser.add_argument('--vmax', type=float, help='Colour scale maximum for mean plots')
    parser.add_argument('--config', help='sweep_config.yaml (default: next to the table)')
    parser.add_argument('--out-dir', help='Figure folder (default: figures/ next to the table)')
    args = parser.parse_args(argv)

    table_dir = os.path.dirname(os.path.abspath(args.table))
    config_path = args.config or os.path.join(table_dir, 'sweep_config.yaml')
    config = load_sweep_config(config_path)

    table = pd.read_csv(args.table)
    out_dir = args.out_dir or os.path.join(table_dir, 'figures')
    for metric in args.metric:
        plot_metric(table, config, metric, out_dir, args.vmin, args.vmax)


if __name__ == '__main__':
    main()
