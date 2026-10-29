"""Plot sweep heatmaps from a results table (made by build_table.py).

    python plot_sweep.py <results.csv> [--metric r_exc r_inh ...] [--axes KEY ...]
                         [--where KEY=VALUE ...] [--vmin V] [--vmax V]
                         [--config PATH] [--out-dir DIR]

For each metric, averages the realizations at every grid point and writes a
mean figure and an SD figure. The plot type depends on how many axes are plotted:
    1 axis  -> line plot with SD error bars
    2 axes  -> heatmap: 1st axis on x, 2nd on y
    4 axes  -> grid of heatmaps: 3rd axis picks the panel row, 4th the panel
               column; inside each panel the 1st axis is x and the 2nd is y

Axis order and labels come from sweep_config.yaml (default: next to the CSV).
Use --axes to plot a subset, and --where to fix the other swept parameters,
e.g. `--axes COUPLING_STRENGTH EPOP_BASE_EXCITE --where G_INTRA=0.5 --where G_INTER=1.0`.
"""
import argparse
import os
import sys
import types
from typing import Dict, List, Sequence, Tuple

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import plotting.sync_heatmap as heatmap
from sweep_config import axis_keys, load_sweep_config

METRIC_LABELS = {
    'r_exc': 'KOP r (E pop)', 'r_inh': 'KOP r (I pop)',
    'chi_exc': r'$\chi$ (E pop)', 'chi_inh': r'$\chi$ (I pop)',
    'rate_exc_hz': 'Firing rate, E pop (Hz)', 'rate_inh_hz': 'Firing rate, I pop (Hz)',
    'median_freq_exc_hz': 'Median frequency, E pop (Hz)',
    'median_freq_inh_hz': 'Median frequency, I pop (Hz)',
    'isi_cv_exc': 'ISI CV (E pop)', 'isi_cv_inh': 'ISI CV (I pop)',
}


def select_rows(table: pd.DataFrame, where: Dict[str, float]) -> pd.DataFrame:
    """Keep rows whose fixed parameters match (to float tolerance)."""
    for key, value in where.items():
        if key not in table.columns:
            sys.exit(f"--where {key}: no such column")
        table = table[np.isclose(table[key], value)]
    if table.empty:
        sys.exit(f"No rows match {where}")
    return table


def check_other_axes_fixed(table: pd.DataFrame, config: dict, plot_keys: List[str]) -> None:
    """Every swept parameter that isn't plotted must have a single value left."""
    for key in axis_keys(config):
        if key in plot_keys:
            continue
        values = sorted(table[key].unique())
        if len(values) > 1:
            sys.exit(f"{key} isn't plotted but has {len(values)} values {values}; "
                     f"fix it with --where {key}=VALUE")


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


def plot_1d(stats, axis, metric_label, out_dir, name):
    mean, (x_values,) = to_grid(stats, 'mean', [axis['key']])
    sd, _ = to_grid(stats, 'std', [axis['key']])
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.errorbar(x_values, mean, yerr=sd, marker='o', capsize=3)
    ax.set_xlabel(axis['label'])
    ax.set_ylabel(metric_label)
    ax.set_title(f'{metric_label}: mean ± SD over realizations', fontweight='bold')
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f'{name}.png')
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved to {path}")


def plot_2d(stats, axes, metric_label, out_dir, name, vmin, vmax, with_sd):
    x, y = axes
    figures = [('mean', metric_label, vmin)]
    if with_sd:
        figures.append(('std', f'SD of {metric_label}', 0))
    for stat, title, lo in figures:
        array, (y_values, x_values) = to_grid(stats, stat, [y['key'], x['key']])
        lo, hi = colour_range(array, lo, vmax if stat == 'mean' else None)
        heatmap.plot_synchrony_single(
            types.SimpleNamespace(figures_dir=out_dir), array, x_values, y_values,
            x['label'], y['label'], title=title, vmin=lo, vmax=hi,
            save_name=f'{name}_{stat}.png', cbar_label=title)


def plot_4d(stats, axes, metric_label, out_dir, name, vmin, vmax, with_sd):
    x, y, row, col = axes
    figures = [('mean', metric_label, vmin)]
    if with_sd:
        figures.append(('std', f'SD of {metric_label}', 0))
    for stat, title, lo in figures:
        array, (row_values, col_values, x_values, y_values) = to_grid(
            stats, stat, [row['key'], col['key'], x['key'], y['key']])
        lo, hi = colour_range(array, lo, vmax if stat == 'mean' else None)
        heatmap.plot_synchrony_multifaceted_single(
            types.SimpleNamespace(figures_dir=out_dir), array,
            row_values, col_values, x_values, y_values,
            row['label'], col['label'], x['label'], y['label'],
            title=title, vmin=lo, vmax=hi,
            save_name=f'{name}_{stat}.png', cbar_label=title)


def plot_metric(table, config, plot_keys, metric, out_dir, vmin, vmax):
    """Write the mean figure, and the SD figure if there are repeat realizations, for one metric."""
    if metric not in table.columns:
        sys.exit(f"--metric {metric}: no such column")
    axes = [next(a for a in config['axes'] if a['key'] == key) for key in plot_keys]
    stats = grid_stats(table, plot_keys, metric)
    label = METRIC_LABELS.get(metric, metric)
    name = f"{metric}_{'_'.join(plot_keys)}"
    # SD needs at least 2 realizations: pandas gives NaN for the SD of one value.
    with_sd = stats['count'].max() > 1
    if not with_sd:
        print(f"  {metric}: one realization per grid point, so no SD figure")
    if len(axes) == 1:
        plot_1d(stats, axes[0], label, out_dir, name)
    elif len(axes) == 2:
        plot_2d(stats, axes, label, out_dir, name, vmin, vmax, with_sd)
    elif len(axes) == 4:
        plot_4d(stats, axes, label, out_dir, name, vmin, vmax, with_sd)
    else:
        sys.exit("Plotting 3 axes isn't supported; plot 2 with --axes and fix the third with --where")
    print(f"  {metric}: jobs per grid point min {stats['count'].min()}, max {stats['count'].max()}")


def parse_where(items: List[str]) -> Dict[str, float]:
    where = {}
    for item in items:
        key, sep, value = item.partition('=')
        if not sep:
            sys.exit(f"--where expects KEY=VALUE, got {item!r}")
        where[key] = float(value)
    return where


def main(argv=None):
    parser = argparse.ArgumentParser(description='Plot sweep heatmaps from results.csv.')
    parser.add_argument('table', help='results.csv from build_table.py')
    parser.add_argument('--metric', nargs='+', default=['r_exc', 'r_inh'],
                        help='Column(s) to plot (default: r_exc r_inh)')
    parser.add_argument('--axes', nargs='+', help='Swept keys to plot, in plot order (default: all)')
    parser.add_argument('--where', action='append', default=[],
                        help='Fix a swept parameter: KEY=VALUE (repeatable)')
    parser.add_argument('--vmin', type=float, help='Colour scale minimum for mean plots')
    parser.add_argument('--vmax', type=float, help='Colour scale maximum for mean plots')
    parser.add_argument('--config', help='sweep_config.yaml (default: next to the table)')
    parser.add_argument('--out-dir', help='Figure folder (default: figures/ next to the table)')
    args = parser.parse_args(argv)

    table_dir = os.path.dirname(os.path.abspath(args.table))
    config_path = args.config or os.path.join(table_dir, 'sweep_config.yaml')
    config = load_sweep_config(config_path)
    plot_keys = args.axes or axis_keys(config)
    unknown = [key for key in plot_keys if key not in axis_keys(config)]
    if unknown:
        sys.exit(f"--axes {unknown} aren't swept in {config_path}")

    table = select_rows(pd.read_csv(args.table), parse_where(args.where))
    check_other_axes_fixed(table, config, plot_keys)
    out_dir = args.out_dir or os.path.join(table_dir, 'figures')
    for metric in args.metric:
        plot_metric(table, config, plot_keys, metric, out_dir, args.vmin, args.vmax)


if __name__ == '__main__':
    main()
