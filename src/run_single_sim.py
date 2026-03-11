#!/usr/bin/env python3

import argparse
import os
import pickle
import numpy as np

# Give each job a unique Brian2 cache directory to avoid C++ compilation
# conflicts between parallel Condor jobs, while keeping the fast C++ backend.
import tempfile
from brian2 import prefs
from brian2 import *

import config
import params
import data_processing
import synch as syn
import plotting as ph
from sim import run_sim


def main():
    parser = argparse.ArgumentParser(description='Run one simulation for a given CE and X0.')
    parser.add_argument('--ce', type=float, required=True, help='Coupling strength CE')
    parser.add_argument('--x0', type=float, required=True, help='Epileptogenicity X0')
    parser.add_argument('--realization', type=int, default=1, help='Realization number (sets random seed)')
    args = parser.parse_args()

    params.COUPLING_STRENGTH = args.ce
    params.HR_X_NAUGHT = args.x0

    job_id = f'CE_{args.ce:.3f}_X0_{args.x0:.3f}_r{args.realization}'

    # Unique cache dir per job prevents parallel Condor jobs from colliding on C++ compilation
    cache_dir = tempfile.mkdtemp(prefix=f'brian2_{job_id}_')
    prefs.codegen.runtime.cython.cache_dir = cache_dir

    # Set random seed for reproducibility across realizations
    seed(args.realization)

    print(f"Starting job: {job_id}")

    # Give this job its own data subdirectory to avoid file conflicts with parallel jobs
    job_data_dir = os.path.join('data', 'jobs', job_id)
    os.makedirs(job_data_dir, exist_ok=True)
    data_processing.DATA_DIR = job_data_dir

    # Run simulation — save_data writes to job_data_dir/output_data.pkl
    run_sim()

    # Load results
    data = data_processing.load_sim_data()
    res = data['results']
    x1 = res['x1']
    x2 = res['x2']
    t  = res['t']

    # Compute synchrony
    chi, _, _ = syn.autocorelate(x1)
    print(f"  chi = {chi:.4f}")

    # Save per-job result to data/results/
    results_dir = os.path.join('data', 'results')
    os.makedirs(results_dir, exist_ok=True)
    job_result = {
        'ce':          args.ce,
        'x0':          args.x0,
        'realization': args.realization,
        'chi':         float(chi),
    }
    with open(os.path.join(results_dir, f'{job_id}.pkl'), 'wb') as f:
        pickle.dump(job_result, f)

    # Save debug plot
    debug_dir = os.path.join('figures', 'sweep_debug')
    os.makedirs(debug_dir, exist_ok=True)
    spike_matrix_1 = data_processing.create_spike_matrix_histo(res['spikes_n1'], params.NUM_CELLS, 0)
    spike_matrix_2 = data_processing.create_spike_matrix_histo(res['spikes_n2'], params.NUM_CELLS, 0)
    ph.standard_plot(t, x1, x2, spike_matrix_1, spike_matrix_2,
                     params.NUM_CELLS, params.SIM_DURATION / second,
                     save_path=os.path.join(debug_dir, f'{job_id}.png'),
                     show=False)

    print(f"Done: {job_id}")


if __name__ == '__main__':
    main()
