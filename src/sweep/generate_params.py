"""Turn a sweep_config.yaml into one params YAML per Condor job.

    python generate_params.py <sweep_config.yaml> --name <run> [--overwrite]

Writes, under src/sweep/params/<run>/:

    sweep_config.yaml        # verbatim copy of the input: what this run sweeps
    jobs/param_1.yaml ...    # one complete params file per job
    job_list.txt             # one job YAML path per line, for Condor's `Queue ... from`

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

from sweep_config import check_keys_in_base, load_sweep_config

SWEEP_DIR = os.path.dirname(os.path.abspath(__file__))
PARAMS_ROOT = os.path.join(SWEEP_DIR, 'params')


def job_grid(config: dict) -> List[Tuple[str, Dict[str, str], int]]:
    """Return a list of (job_id, swept_values, seed), one entry per job.

    swept_values holds this job's value for each swept parameter, e.g.
    {'COUPLING_STRENGTH': '0.05*1', 'G_INTRA': '0.5*uS'}.

    Jobs are numbered from 1. The first axis changes slowest and the
    realization changes fastest.
    """

    # get values for each sweep axis into a list of lists
    axis_values = [axis['values'] for axis in config['axes']]
    seeds = range(1, config['realizations'] + 1)

    jobs = []
    job_number = 1
    # Each grid point is a tuple with one value per axis, e.g. (0.05, 0.5).
    # points are computed by unpacking the axis value lists and taking the cartesian product of their elements
    for grid_point in product(*axis_values):
        swept_values = {}
        # runs for the number of dimensions swept over.
        # eg: for a 4d plot this runs 4 times
        for axis, value in zip(config['axes'], grid_point):
            # forms a dict of the values being swept over
            # grid_point is purely a tuple of values
            # add back the key and unit to the value
            swept_values[axis['key']] = format_value(value, axis['unit'])

        # One job per realization at this grid point.
        for seed in seeds:
            jobs.append((f'param_{job_number}', swept_values, seed))
            job_number += 1
    return jobs


def format_value(value: float, unit: str) -> str:
    """Write a swept value the way params YAMLs do, e.g. '0.5*uS' or '-4.5*1'."""
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
    """Create an empty run directory, refusing to reuse a name unless asked."""
    if os.path.exists(run_dir):
        if not overwrite:
            sys.exit(f"Run directory {run_dir} already exists. Pick a new --name, "
                     f"or pass --overwrite (only if none of its jobs are still queued).")
        shutil.rmtree(run_dir)
    os.makedirs(os.path.join(run_dir, 'jobs'))


def write_jobs(config: dict, run_dir: str) -> List[str]:
    """Write every job YAML; return their paths relative to SWEEP_DIR."""
    with open(config['base_params']) as f:
        base = yaml.safe_load(f)

    job_paths = []
    for job_id, swept_values, seed in job_grid(config):
        path = os.path.join(run_dir, 'jobs', f'{job_id}.yaml')
        with open(path, 'w') as f:
            yaml.safe_dump(job_params(base, swept_values, seed, config), f, sort_keys=False)
        # Relative to src/sweep/, which is the Condor Initialdir and the
        # directory wrapper.sh runs jobs from.
        job_paths.append(os.path.relpath(path, SWEEP_DIR))
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
    parser.add_argument('--name', '-n', required=True, help='Run name (folder under params/)')
    parser.add_argument('--overwrite', action='store_true',
                        help='Replace an existing run with the same name')
    args = parser.parse_args(argv)

    config = load_sweep_config(args.config)
    check_keys_in_base(config)

    run_dir = os.path.join(PARAMS_ROOT, args.name)
    prepare_run_dir(run_dir, args.overwrite)
    shutil.copyfile(args.config, os.path.join(run_dir, 'sweep_config.yaml'))
    job_paths = write_jobs(config, run_dir)
    list_path = write_job_list(job_paths, run_dir)

    shape = ' x '.join(f"{len(axis['values'])} {axis['key']}" for axis in config['axes'])
    print(f"Grid: {shape} x {config['realizations']} realizations = {len(job_paths)} jobs")
    print(f"Wrote {run_dir}/jobs/ and {list_path}")


if __name__ == '__main__':
    main()
