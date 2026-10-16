#!/usr/bin/env python3
"""Run one sweep job: simulate one job YAML and save its raw output.

    python run_single_sim.py --params params/<run>/jobs/param_7.yaml [--out-dir DIR] [--debug-plot]

Condor runs one copy per line of job_list.txt. No metrics are computed here;
build_table.py derives them from the raw files, so they can be recomputed later.

Output (same layout as older full-output pickles, minus the synaptic current traces):
    {'metadata': {..., job_id, cb_on, low_res, seed},
     'params':   the job's params dict,
     'results':  {t, x1, x2, spikes_n1: {t, i}, spikes_n2: {t, i}}}

Where it goes:
    - Under Condor: $_CONDOR_SCRATCH_DIR/output_<job_id>.pkl. Condor copies it to
      output/sweep/<run>/jobs/<job_id>/output.pkl (see condor/setup_condor.sh).
    - Locally: <out-dir>/jobs/<job_id>/output.pkl.
"""
import argparse
import os
import pickle
import sys
import tempfile
import types

from brian2 import prefs, seed, second

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from model import run_sim
from param_loader import load_params

RAW_RESULT_KEYS = ('t', 'x1', 'x2', 'spikes_n1', 'spikes_n2')


def job_settings(params: dict, params_path: str) -> dict:
    """Read the per-job settings that generate_params.py writes into every job YAML."""
    missing = [k for k in ('SEED', 'CB_ON', 'LOW_RES') if k not in params]
    if missing:
        sys.exit(f"{params_path} has no {missing}; make job YAMLs with generate_params.py")
    return {
        'job_id': os.path.splitext(os.path.basename(params_path))[0],
        'seed': int(params['SEED']),
        'cb_on': bool(params['CB_ON']),
        'low_res': bool(params['LOW_RES']),
    }


def simulate(params: dict, settings: dict) -> dict:
    """Run the model in memory and return its full sim_data dict."""
    # Separate Cython cache per job so parallel jobs don't collide while compiling.
    prefs.codegen.runtime.cython.cache_dir = tempfile.mkdtemp(prefix=f"brian2_{settings['job_id']}_")
    seed(settings['seed'])
    return run_sim(None, params, cb_on=settings['cb_on'], save=False,
                   low_res=settings['low_res'])


def raw_output(sim_data: dict, settings: dict) -> dict:
    """Keep only what metrics need, and record how the job was run."""
    return {
        'metadata': {**sim_data['metadata'], **settings},
        'params': sim_data['params'],
        'results': {key: sim_data['results'][key] for key in RAW_RESULT_KEYS},
    }


def output_path(job_id: str, out_dir: str) -> str:
    """Condor scratch root if running under Condor, else <out_dir>/jobs/<job_id>/output.pkl."""
    scratch = os.environ.get('_CONDOR_SCRATCH_DIR')
    if scratch:
        return os.path.join(scratch, f'output_{job_id}.pkl')
    if out_dir is None:
        sys.exit("Not running under Condor: pass --out-dir for local runs")
    job_dir = os.path.join(out_dir, 'jobs', job_id)
    os.makedirs(job_dir, exist_ok=True)
    return os.path.join(job_dir, 'output.pkl')


def save_debug_plot(sim_data: dict, params: dict, figures_dir: str) -> None:
    """Write the standard population plot. Local runs only: Condor doesn't copy it back."""
    import matplotlib
    matplotlib.use('Agg')
    import numpy as np
    import data_processing
    import plotting.population_plots as population_plots

    os.makedirs(figures_dir, exist_ok=True)
    res = sim_data['results']
    num_cells = params['NUM_CELLS']
    bin_width = params['PLOT_BIN_WIDTH'] / second
    histo_bin_width = params['SPIKE_HISTO_BIN_WIDTH'] / second
    population_plots.standard_plot(
        types.SimpleNamespace(data_dir=None, figures_dir=figures_dir), params,
        res['t'],
        data_processing.bin_data(res['t'], np.mean(res['x1'], axis=0), bin_width),
        data_processing.bin_data(res['t'], np.mean(res['x2'], axis=0), bin_width),
        data_processing.create_spike_matrix_histo(params, res['spikes_n1'], num_cells, histo_bin_width),
        data_processing.create_spike_matrix_histo(params, res['spikes_n2'], num_cells, histo_bin_width),
        num_cells, params['SIM_DURATION'] / second,
        g_inter_vals=params['G_INTER_VALS'], g_intra_vals=params['G_INTRA_VALS'],
        x0_t=res.get('x0_t'), ce_t=res.get('ce_t'), show=False)


def main(argv=None):
    parser = argparse.ArgumentParser(description='Run one sweep job and save its raw output.')
    parser.add_argument('--params', required=True,
                        help='Job YAML written by generate_params.py')
    parser.add_argument('--out-dir', default=None,
                        help='Local runs only: write <out-dir>/jobs/<job_id>/output.pkl')
    parser.add_argument('--debug-plot', action='store_true',
                        help='Local runs only: also write <out-dir>/figures/<job_id>/ plots')
    args = parser.parse_args(argv)

    params = load_params(args.params)
    settings = job_settings(params, args.params)
    path = output_path(settings['job_id'], args.out_dir)

    print(f"Starting job {settings['job_id']} (seed {settings['seed']}, cb_on {settings['cb_on']})")
    sim_data = simulate(params, settings)
    with open(path, 'wb') as f:
        pickle.dump(raw_output(sim_data, settings), f)
    print(f"Wrote {path}")

    if args.debug_plot and args.out_dir is not None:
        save_debug_plot(sim_data, params, os.path.join(args.out_dir, 'figures', settings['job_id']))


if __name__ == '__main__':
    main()
