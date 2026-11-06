"""Build a results table (one row per job) from a sweep's raw output files.

    python build_table.py <run_dir> [--config PATH] [--out PATH]

Reads <run_dir>/runs/<job_id>.pkl, computes every summary metric with synch_metrics.summarize_run, and
writes a CSV:

    job_id, <one column per swept key, in the config's units>, seed, cb_on, ...,
    chi_exc, chi_inh, r_exc, r_inh, rate_exc_hz, ..., source_path

The raw files are never modified, so the table can be rebuilt at any time,
e.g. after a metric fix.
"""
import argparse
import glob
import os
import pickle
import sys
from typing import List

import brian2
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import synch_metrics as syn
from sweep_config import load_sweep_config, sweep_param_keys

def find_raw_files(run_dir: str) -> List[str]:
    """All runs/<job_id>.pkl files under run_dir, sorted by path."""
    return sorted(glob.glob(os.path.join(run_dir, 'runs', '*.pkl')))


def to_unit(value, unit: str) -> float:
    """Express a value in the declared unit, e.g. 5e-7 S -> 0.5 for unit 'uS'."""
    # The unit is a string like 'uS' or '1'; look it up among Brian's names.
    ratio = value / eval(str(unit), {}, vars(brian2))
    if isinstance(ratio, brian2.Quantity) and not ratio.is_dimensionless:
        raise ValueError(f"{value} can't be expressed in {unit}")
    number = float(ratio)
    # Strip float noise from the unit division so equal values group together.
    return float(f'{number:.12g}')

# why is this needed
def to_seconds(value) -> float:
    """Seconds as a float, whether stored as a Brian quantity or a bare number."""
    if isinstance(value, brian2.Quantity):
        return float(value / brian2.second)
    return float(value)


def create_job_row(path: str, config: dict) -> dict:
    """One table row: the job's swept values, run settings, and metrics."""
    with open(path, 'rb') as f:
        raw = pickle.load(f)
    params, results = raw['params'], raw['results']
    metadata = raw['metadata']

    duration = to_seconds(params['SIM_DURATION'])
    t_start = to_seconds(params['TRANSIENT'])
    num_cells = int(params['NUM_CELLS'])

    row = {'job_id': metadata['job_id']}
    for sweep_param in config['sweep_params']:
        key = sweep_param['key']
        row[key] = to_unit(params[key], sweep_param['unit'])
    row.update({
        'seed': metadata['seed'],
        'cb_on': metadata['cb_on'],
        'low_res': metadata['low_res'],
        'num_cells': num_cells,
        'sim_duration_s': duration,
    })

    row.update(syn.summarize_run(results, num_cells, duration, t_start))
    row['r_sum'] = row['r_inh'] + row['r_exc']
    row['source_path'] = os.path.abspath(path)
    return row


def build_table(run_dir: str, config: dict) -> pd.DataFrame:
    """Rows for every raw file in run_dir. Files that fail to load are reported and skipped."""
    paths = find_raw_files(run_dir)
    if not paths:
        sys.exit(f"No runs/*.pkl files under {run_dir}")

    rows, failed = [], []
    for n, path in enumerate(paths, start=1):
        print(f"[{n}/{len(paths)}] {path}")
        try:
            rows.append(create_job_row(path, config))
        except Exception as err:  # keep going; report at the end
            failed.append((path, f'{type(err).__name__}: {err}'))
            
    for path, reason in failed:
        print(f"FAILED {path}: {reason}")
    if not rows:
        sys.exit("Every file failed; no table written")
    return pd.DataFrame(rows).sort_values(sweep_param_keys(config) + ['seed'], ignore_index=True)


def print_grid_summary(table: pd.DataFrame, config: dict) -> None:
    """Print each swept parameter's values, so a wrong config or a broken sweep is obvious."""
    print(f"\n{len(table)} rows")
    for key in sweep_param_keys(config):
        values = sorted(float(v) for v in table[key].unique())
        print(f"  {key}: {len(values)} values {values}")
    counts = table.groupby(sweep_param_keys(config)).size()
    print(f"  jobs per grid point: min {counts.min()}, max {counts.max()}")


def main(argv=None):
    parser = argparse.ArgumentParser(description='Build results.csv from a sweep\'s raw output.')
    parser.add_argument('run_dir', help='Sweep folder, data/<sweep>/')
    parser.add_argument('--config', help='sweep_config.yaml (default: <run_dir>/sweep_config.yaml)')
    parser.add_argument('--out', help='Output CSV (default: <run_dir>/results.csv)')
    args = parser.parse_args(argv)

    config_path = args.config or os.path.join(args.run_dir, 'sweep_config.yaml')
    try:
        config = load_sweep_config(config_path)
    except (OSError, ValueError) as err:
        sys.exit(f"Can't read sweep config {config_path}: {err}")

    table = build_table(args.run_dir, config)
    out = args.out or os.path.join(args.run_dir, 'results.csv')
    table.to_csv(out, index=False)
    print_grid_summary(table, config)
    print(f"Wrote {out}")


if __name__ == '__main__':
    main()
