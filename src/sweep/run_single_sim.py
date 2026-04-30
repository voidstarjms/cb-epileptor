#!/usr/bin/env python3
"""One sim job in the parameter sweep — Condor runs one process per YAML in params/.

Each job's full configuration (including sweep dimensions and seed) lives in the
YAML pointed to by --params. Outputs go under cwd/data and cwd/figures, so this
should be run with cwd = src/sweep/ (Condor's Initialdir handles that).
"""

import argparse
import os
import pickle
import tempfile

# Headless backend. Must be set before any pyplot import, direct or transitive.
# import matplotlib
# matplotlib.use('Agg')

from brian2 import prefs, seed, second, uS

import sys
sys.path.append("..")
import data_processing
import synch as syn
import plotting.population_plots as ph
from model import run_sim
from param_loader import load_params
from run import FilePaths


def _to_uS_float(v):
    """Convert a Brian2 conductance quantity (or plain number) to a float in uS.

    Args:
        v: Brian2 conductance quantity (e.g. 0.5 * uS) or a plain number.

    Returns:
        float: The value expressed in microsiemens.
    """
    try:
        return float(v / uS)
    except (TypeError, ValueError):
        return float(v)


def main():
    """Run one sim from a YAML, write its chi summary plus a debug plot.

    Reads the YAML pointed to by --params, runs the simulation via model.run_sim,
    computes synchrony chi over the HR population, and writes:

      - data/jobs/<job_id>/output.pkl: full sim output.
      - data/results/<job_id>.pkl: compact summary {ce, x0, Gintra, Ginter,
        realization, chi} consumed by aggregate.py.
      - figures/sweep_debug/<job_id>/standard_plot.png: per-job debug plot.
    """
    parser = argparse.ArgumentParser(description='Run one simulation from a YAML param file.')
    parser.add_argument('--params', type=str, required=True,
                        help='Path to a YAML param file (e.g. params/param_1.yaml).')
    args = parser.parse_args()

    params_dict = load_params(args.params)
    job_id = os.path.splitext(os.path.basename(args.params))[0]

    # Unique cache dir so parallel jobs don't collide on C++ compilation.
    cache_dir = tempfile.mkdtemp(prefix=f'brian2_{job_id}_')
    prefs.codegen.runtime.cython.cache_dir = cache_dir

    realization = int(params_dict.get('SEED', 1))
    seed(realization)

    print(f"Starting job: {job_id}")

    # Per-job dirs so parallel jobs don't overwrite each other's output.
    filepaths = FilePaths(
        data_dir=os.path.join('data', 'jobs', job_id),
        figures_dir=os.path.join('figures', 'sweep_debug', job_id),
    )
    os.makedirs(filepaths.data_dir, exist_ok=True)
    os.makedirs(filepaths.figures_dir, exist_ok=True)

    run_sim(filepaths, params_dict)

    data = data_processing.load_sim_data(filepaths)
    res = data['results']
    x1 = res['x1']
    x2 = res['x2']
    t  = res['t']

    # Compute synchrony
    chi, _, _ = syn.autocorrelate(x1)
    print(f"  chi = {chi:.4f}")

    # Pull the sweep dims back out of the YAML so aggregate.py can reconstruct
    # the grid without re-parsing filenames.
    ce     = float(params_dict['COUPLING_STRENGTH'])
    x0     = float(params_dict['HR_X_NAUGHT'])
    Gintra = _to_uS_float(params_dict['G_INTRA'])
    Ginter = _to_uS_float(params_dict['G_INTER'])

    # Save per-job result to data/results/
    results_dir = os.path.join('data', 'results')
    os.makedirs(results_dir, exist_ok=True)
    job_result = {
        'ce':          ce,
        'x0':          x0,
        'Gintra':      Gintra,
        'Ginter':      Ginter,
        'realization': realization,
        'chi':         float(chi),
    }
    with open(os.path.join(results_dir, f'{job_id}.pkl'), 'wb') as f:
        pickle.dump(job_result, f)

    spike_matrix_1 = data_processing.create_spike_matrix_histo(
        params_dict, res['spikes_n1'], params_dict['NUM_CELLS'])
    spike_matrix_2 = data_processing.create_spike_matrix_histo(
        params_dict, res['spikes_n2'], params_dict['NUM_CELLS'])
    ph.standard_plot(filepaths, params_dict, t, x1, x2,
                     spike_matrix_1, spike_matrix_2,
                     params_dict['NUM_CELLS'],
                     params_dict['SIM_DURATION'] / second,
                     g_inter_vals=params_dict['G_INTER_VALS'],
                     g_intra_vals=params_dict['G_INTRA_VALS'],
                     coupling_vals=params_dict['COUPLING_VALS'],
                     x_naught_vals=params_dict['X_NAUGHT_VALS'],
                     show=False)

    print(f"Done: {job_id}")


if __name__ == '__main__':
    main()
