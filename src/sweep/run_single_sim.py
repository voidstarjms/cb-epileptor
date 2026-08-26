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
import numpy as np

# Headless backend. Must be set before any pyplot import, direct or transitive.
# import matplotlib
# matplotlib.use('Agg')

from brian2 import prefs, seed
from brian2.units import *

import sys
sys.path.append("..")
import data_processing
import synch_metrics as syn
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
    """Run one sim from a YAML, write its chi summary; optional per-job debug plot.

    Reads the YAML pointed to by --params and , runs the simulation via model.run_sim
    in memory (no intermediate output.pkl on NFS), computes synchrony chi and
    KOP r over the HR population, and writes:

      - data/results/metrics_<job_id>.pkl: compact summary {ce, x0, Gintra, Ginter,
        realization, chi, r} consumed by aggregate.py. (Only NFS write per job.)
      - figures/sweep_debug/<job_id>/standard_plot.svg: per-job debug plot,
        written only when --debug-plot is passed.
    """
    parser = argparse.ArgumentParser(description='Run one simulation from a YAML param file.')
    parser.add_argument('--params', type=str, default='../../params.yaml', required=False,
                        help='Path to a YAML param file (e.g. params/param_1.yaml).')
    parser.add_argument('--excite', type=float, default=None, required=False,
                        help='Manually assigns epop excitability')
    parser.add_argument('--J', type=float, default=None, required=False,
                        help='Manually assigns gap junction strength')
    parser.add_argument('--Gintra', type=float, default=None, required=False,
                        help='Manually assigns intrapopulation synapse strength')
    parser.add_argument('--Ginter', type=float, default=None, required=False,
                        help='Manually assigns interpopulation synapse strength')
    parser.add_argument('--realization', type=int, default=None, required=False,
                        help='The random seed for the experiment')
    parser.add_argument('--debug-plot', action='store_true', default=False,
                        help='Write a per-job standard_plot debug image. Off by default '
                             'for sweeps to avoid per-job NFS plot writes.')
    args = parser.parse_args()
    
    params_dict = load_params(args.params)
    # Load scientific parameters from command line if supplied
    if args.excite != None:
        params_dict['EPOP_BASE_EXCITE'] = args.excite
    if args.J != None:
        params_dict['COUPLING_STRENGTH'] = args.J
    if args.Gintra != None:
        params_dict['G_INTRA'] = args.Gintra
    if args.Ginter != None:
        params_dict['G_INTER'] = args.Ginter
    if args.realization != None:
        params_dict['SEED'] = args.realization

    # Job ID is the YAML basename (e.g. 'param_37'). Matches what condor/README.md documents.
    job_id = os.path.splitext(os.path.basename(args.params))[0]

    # Unique cache dir so parallel jobs don't collide on C++ compilation.
    cache_dir = tempfile.mkdtemp(prefix=f'brian2_{job_id}_')
    prefs.codegen.runtime.cython.cache_dir = cache_dir

    realization = int(params_dict.get('SEED', args.realization))
    seed(realization)

    print(f"Starting job: {job_id}")

    # Per-job dirs so parallel jobs don't overwrite each other's output.
    # data_dir is unused under the in-memory pipeline but kept on the dataclass
    # for plot helpers that expect a populated FilePaths.
    filepaths = FilePaths(
        data_dir=os.path.join('data', 'jobs', job_id),
        figures_dir=os.path.join('figures', 'sweep_debug', job_id),
    )

    # In-memory pipeline: no output.pkl write/read/delete on NFS.
    data = run_sim(filepaths, params_dict, cb_on=True, save=False)
    res = data['results']
    x1 = res['x1']
    x2 = res['x2']
    t  = res['t']
    spikes_1 = res['spikes_n1']
    spikes_2 = res['spikes_n2']

    num_cells = params_dict['NUM_CELLS']
    spike_matrix_1 = data_processing.create_spike_matrix_histo(params_dict, res['spikes_n1'], num_cells)
    spike_matrix_2 = data_processing.create_spike_matrix_histo(params_dict, res['spikes_n2'], num_cells)

    # Compute synchrony
    chi_exc, _, _ = syn.autocorrelate(x1)
    chi_inh, _, _ = syn.autocorrelate(x2)

    _, r_exc, _ = syn.KOP(spikes_1['i'], spikes_1['t'], params_dict['SIM_DURATION'] / second)
    _, r_inh, _ = syn.KOP(spikes_2['i'], spikes_2['t'], params_dict['SIM_DURATION'] / second)
    r = np.mean(r_exc) + np.mean(r_inh)
    
    print(f"  r = {r:.4f}")

    # Pull the sweep dims back out of the YAML so aggregate.py can reconstruct
    # the grid without re-parsing filenames.
    coupling_J = float(params_dict['COUPLING_STRENGTH'])
    epop_excite = float(params_dict['EPOP_BASE_EXCITE'])
    Gintra = _to_uS_float(params_dict['G_INTRA'])
    Ginter = _to_uS_float(params_dict['G_INTER'])

    # Save per-job result. Under Condor file transfer the job runs in a scratch
    # sandbox; write to its root ($_CONDOR_SCRATCH_DIR) so Condor transfers the metrics
    # file back. Falls back to data/results/ for local (non-Condor) runs.
    results_dir = os.environ.get('_CONDOR_SCRATCH_DIR', os.path.join('data', 'results'))
    os.makedirs(results_dir, exist_ok=True)
    job_result = {
        'params':           params_dict,
        'realization':      realization,
        'chi_exc':          float(chi_exc),
        'chi_inh':          float(chi_inh),
        'r_exc':            float(np.mean(r_exc)),
        'r_inh':            float(np.mean(r_inh)),
        'exc_spike_mat':    spike_matrix_1.astype(np.uint8),
        'inh_spike_mat':    spike_matrix_2.astype(np.uint8),
    }
    
    with open(os.path.join(results_dir, f'metrics_{job_id}.pkl'), 'wb') as f:
        pickle.dump(job_result, f)

    if args.debug_plot:
        os.makedirs(filepaths.figures_dir, exist_ok=True)
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
                         x0_t=res.get('x0_t'), ce_t=res.get('ce_t'),
                         show=False)

    print(f"Done: {job_id}")


if __name__ == '__main__':
    main()
