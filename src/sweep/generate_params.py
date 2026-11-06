"""Turn a sweep_config.yaml into one params YAML per Condor job.

    python generate_params.py <sweep_config.yaml> --name <run> [--overwrite]

Writes, under data/sweep/<run>/:

    sweep_config.yaml        # verbatim copy of the input: what this run sweeps
    base_params.yaml         # verbatim copy of the base params every job starts from
    jobs/param_1.yaml ...    # one complete params file per job
    job_list.txt             # one job YAML path per line, for Condor's `Queue ... from`

The same folder later holds runs/ (one raw .pkl per job), results.csv and figures/.

Each job YAML is the base params file with the swept keys replaced (e.g.
G_INTRA: "0.5*uS"), plus SEED, CB_ON and LOW_RES. One job per
(grid point, realization); realization k uses SEED = k at every grid point.
"""
import argparse
from itertools import product
import os
import shutil
import sys
from typing import Dict, List, Tuple

import yaml

from sweep_config import load_sweep_config, validate_sweep_params

PROJ_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SWEEP_DATA_ROOT = os.path.join(PROJ_DIR, 'data', 'sweep')


def job_grid(config: dict) -> List[Tuple[str, Dict[str, str], int]]:
    """Return a list of (job_id, swept_values, seed), one entry per job.

    swept_values holds this job's value for each swept parameter, e.g.
    {'COUPLING_STRENGTH': '0.05*1', 'G_INTRA': '0.5*uS'}.

    Jobs are numbered from 1. The first swept parameter changes slowest and the
    realization changes fastest.
    """

    # get values for each swept parameter into a list of lists
    sweep_values = [p['values'] for p in config['sweep_params']]

    jobs = []
    job_number = 1 # counter
    # Each grid point is a tuple with one value per swept parameter, e.g. (0.05, 0.5).
    # points are computed by unpacking the value lists and taking the cartesian product of their elements
    for grid_point in product(*sweep_values):
        # grid point is a tuple of the sweep values. eg (ce, x0, gintra, ginter)
        swept_values = {}
        # runs for the number of dimensions swept over.
        # eg: for a 4d plot this runs 4 times
        for sweep_param, value in zip(config['sweep_params'], grid_point):
            # forms a dict of the values being swept over
            # grid_point is purely a tuple of values
            # add back the key and unit to the value
            swept_values[sweep_param['key']] = format_value(value, sweep_param['unit'])

        # One job per realization at this grid point; realization k uses SEED = k.
        for seed in range(1, config['realizations'] + 1):
            jobs.append((f'param_{job_number}', swept_values, seed))
            job_number += 1
    return jobs


def format_value(value: float, unit: str) -> str:
    """Attach brian units to a value"""
    return f'{float(value)}*{unit}'


def job_params(base: dict, swept_values: Dict[str, str], seed: int, config: dict) -> dict:
    """Return a copy of the base params with this job's values filled in."""
    params = dict(base)
    params.update(swept_values)
    params['SEED'] = seed
    params['CB_ON'] = config['cb_on']
    params['LOW_RES'] = config['low_res']
    return params


def prepare_run_dir(run_dir: str, overwrite: bool) -> None:
    """Create an empty run directory."""
    if os.path.exists(run_dir):
        if not overwrite:
            sys.exit(f"Run directory {run_dir} already exists. Pick a new --name, "
                     f"or pass --overwrite (only if none of its jobs are still queued).")
        # The run directory also holds the simulation results; never delete those.
        runs_dir = os.path.join(run_dir, 'runs')
        if os.path.isdir(runs_dir) and os.listdir(runs_dir):
            sys.exit(f"{runs_dir} already has results; --overwrite won't delete them. "
                     f"Pick a new --name, or remove the folder yourself.")
        shutil.rmtree(run_dir)
    os.makedirs(os.path.join(run_dir, 'jobs'))


def write_jobs(config: dict, run_dir: str) -> List[str]:
    """Write every job YAML; return their paths relative to PROJ_DIR."""
    with open(config['base_params']) as f:
        base = yaml.safe_load(f)

    job_paths = []
    for job_id, swept_values, seed in job_grid(config):
        path = os.path.join(run_dir, 'jobs', f'{job_id}.yaml')
        with open(path, 'w') as f:
            yaml.safe_dump(job_params(base, swept_values, seed, config), f, sort_keys=False)
        # Relative to the project root, which is the Condor Initialdir.
        job_paths.append(os.path.relpath(path, PROJ_DIR))
    return job_paths


def write_job_list(job_paths: List[str], run_dir: str) -> str:
    """Write job_list.txt (one job YAML path per line); return its path."""
    path = os.path.join(run_dir, 'job_list.txt')
    with open(path, 'w') as f:
        f.write('\n'.join(job_paths) + '\n')
    return path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('config', help='Path to the sweep_config.yaml')
    parser.add_argument('--name', '-n', required=True, help='Run name (folder under data/sweep/)')
    parser.add_argument('--overwrite', action='store_true',
                        help='Replace an existing run with the same name')
    args = parser.parse_args(argv)

    config = load_sweep_config(args.config)
    validate_sweep_params(config)

    run_dir = os.path.join(SWEEP_DATA_ROOT, args.name)
    prepare_run_dir(run_dir, args.overwrite)
    shutil.copyfile(args.config, os.path.join(run_dir, 'sweep_config.yaml'))
    shutil.copyfile(config['base_params'], os.path.join(run_dir, 'base_params.yaml'))
    job_paths = write_jobs(config, run_dir)
    list_path = write_job_list(job_paths, run_dir)

    shape = ' x '.join(f"{len(p['values'])} {p['key']}" for p in config['sweep_params'])
    print(f"Grid: {shape} x {config['realizations']} realizations = {len(job_paths)} jobs")
    print(f"Wrote {run_dir}/jobs/ and {list_path}")


if __name__ == '__main__':
    main()
