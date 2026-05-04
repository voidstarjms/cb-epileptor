#!/usr/bin/env python3
"""One sim job in the CE x X0 sweep. Run by condor, one process per (ce, x0, realization)."""

import argparse
import os
import pickle
import tempfile
import numpy as np

# Headless backend. Must be set before any pyplot import, direct or transitive.
# import matplotlib
# matplotlib.use('Agg')

from brian2 import prefs, seed
from brian2.units import *

import sys
sys.path.append("..")
import data_processing
import synch as syn
import plotting.population_plots as ph
from model import run_sim
from param_loader import load_params
from run import FilePaths


def main():
    """Run one sim and write its chi summary plus a debug plot."""
    parser = argparse.ArgumentParser(description='Run one simulation for a given CE and X0.')
    parser.add_argument('--ce', type=float, required=True, help='Coupling strength CE')
    parser.add_argument('--Gintra', type=float, required=True, help='Intrapopulation synapse strength')
    parser.add_argument('--Ginter', type=float, required=True, help='Interpopulation synapse strength')
    parser.add_argument('--x0', type=float, required=True, help='Epileptogenicity X0')
    parser.add_argument('--realization', type=int, default=1,
                        help='Realization number (sets random seed)')
    parser.add_argument('--params', type=str, default='../../params.yaml',
                        help="YAML params file (default: '../../params.yaml')")
    args = parser.parse_args()

    params_dict = load_params(args.params)
    # Collapse the sweep dimensions to one value each for this job.
    params_dict['COUPLING_STRENGTH'] = args.ce
    params_dict['HR_X_NAUGHT'] = args.x0
    params_dict['G_INTRA'] = args.Gintra * usiemens
    params_dict['G_INTER'] = args.Ginter * usiemens
    
    job_id = f'CE_{args.ce:.3f}_X0_{args.x0:.3f}_inter_{args.Ginter}_intra_{args.Gintra}_r{args.realization}'

    # Unique cache dir so parallel jobs don't collide on C++ compilation.
    cache_dir = tempfile.mkdtemp(prefix=f'brian2_{job_id}_')
    prefs.codegen.runtime.cython.cache_dir = cache_dir

    seed(args.realization)

    print(f"Starting job: {job_id}")

    # Per-job dirs so parallel jobs don't overwrite each other's output.
    filepaths = FilePaths(
        data_dir=os.path.join('data', 'jobs', job_id),
        figures_dir=os.path.join('figures', 'sweep_debug', job_id),
    )
    os.makedirs(filepaths.data_dir, exist_ok=True)
    os.makedirs(filepaths.figures_dir, exist_ok=True)

    run_sim(filepaths, params_dict, cb_on=False)

    data = data_processing.load_sim_data(filepaths)
    res = data['results']
    x1 = res['x1']
    x2 = res['x2']
    t  = res['t']
    spikes_1 = res['spikes_n1']

    # Compute synchrony
    chi, _, _ = syn.autocorrelate(x1)
    _, r, _ = syn.KOP(spikes_1['i'], spikes_1['t'], params_dict['SIM_DURATION'] / second)
    r = np.mean(r)
    print(f"  chi = {chi:.4f}")
    print(f"  r = {r:.4f}")

    # Save per-job result to data/results/
    results_dir = os.path.join('data', 'results')
    os.makedirs(results_dir, exist_ok=True)
    job_result = {
        'ce':          args.ce,
        'x0':          args.x0,
        'Gintra':      args.Gintra,
        'Ginter':      args.Ginter,
        'realization': args.realization,
        'chi':         float(chi),
        'r':           float(r),
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
