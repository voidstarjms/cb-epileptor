"""Build a results table (one row per job) from a sweep's raw output files.

    python build_table.py <run_output_dir> [--config PATH] [--out PATH] [--cb-on true|false]

Reads <run_output_dir>/jobs/<job_id>/output.pkl (new and old sweeps use this
layout), computes every summary metric with synch_metrics.summarize_run, and
writes a CSV:

    job_id, <one column per swept key, in the config's units>, seed, cb_on, ...,
    chi_exc, chi_inh, r_exc, r_inh, rate_exc_hz, ..., source_path

The raw files are never modified, so the table can be rebuilt at any time,
e.g. after a metric fix. For an old sweep, write a sweep_config.yaml by hand
(only key, unit and label are needed per axis) and pass it with --config.
"""
import argparse
import glob
import math
import os
import pickle
import sys
import warnings
from typing import List, Optional

import brian2
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import synch_metrics as syn
from sweep_config import axis_keys, load_sweep_config

# For these keys the model reads a schedule list, not the scalar key, so the
# schedule holds the value the simulation actually used. Older runs named the
# x0 schedule X_NAUGHT_VALS.
SCHEDULE_KEYS = {
    'EPOP_BASE_EXCITE': ('BASE_EXCITE_VALS', 'X_NAUGHT_VALS'),
    'HR_X_NAUGHT': ('X_NAUGHT_VALS', 'BASE_EXCITE_VALS'),
    'COUPLING_STRENGTH': ('COUPLING_VALS',),
    'G_INTRA': ('G_INTRA_VALS',),
    'G_INTER': ('G_INTER_VALS',),
}


def find_raw_files(run_dir: str) -> List[str]:
    """All jobs/<job_id>/output.pkl files under run_dir, sorted by path."""
    return sorted(glob.glob(os.path.join(run_dir, 'jobs', '*', 'output.pkl')))


def used_value(params: dict, key: str):
    """The value of `key` the simulation actually used.

    Returns the schedule's single entry if the key has one, else the scalar.
    Returns NaN (with a warning) if the schedule has more than one entry: that
    job didn't hold the parameter at one value, so it has no place on the axis.
    """
    for schedule_key in SCHEDULE_KEYS.get(key, ()):
        if schedule_key in params:
            schedule = params[schedule_key]
            if len(schedule) != 1:
                warnings.warn(f"{schedule_key} has {len(schedule)} entries ({schedule}); "
                              f"{key} was not held at one value in this job")
                return math.nan
            return schedule[0]
    if key not in params:
        raise KeyError(f"{key} not found in the job's params")
    return params[key]


def in_axis_unit(value, axis: dict) -> float:
    """Express a value in the axis's declared unit, e.g. 5e-7 S -> 0.5 for unit 'uS'."""
    if isinstance(value, float) and math.isnan(value):
        return value
    # The unit is a string like 'uS' or '1'; look it up among Brian's names.
    unit = eval(str(axis['unit']), {}, vars(brian2))
    ratio = value / unit
    if isinstance(ratio, brian2.Quantity) and not ratio.is_dimensionless:
        raise ValueError(f"{axis['key']} = {value} can't be expressed in {axis['unit']}")
    number = float(ratio)
    # Strip float noise from the unit division so equal values group together.
    return float(f'{number:.12g}')


def to_seconds(value) -> float:
    """Seconds as a float, whether stored as a Brian quantity or a bare number."""
    if isinstance(value, brian2.Quantity):
        return float(value / brian2.second)
    return float(value)


def job_row(path: str, config: dict, cb_on_override: Optional[bool]) -> dict:
    """One table row: the job's swept values, run settings, and metrics."""
    with open(path, 'rb') as f:
        raw = pickle.load(f)
    params, results = raw['params'], raw['results']
    metadata = raw.get('metadata', {})

    duration = to_seconds(params['SIM_DURATION'])
    t_start = to_seconds(params['TRANSIENT'])
    num_cells = int(params['NUM_CELLS'])

    row = {'job_id': os.path.basename(os.path.dirname(path))}
    for axis in config['axes']:
        row[axis['key']] = in_axis_unit(used_value(params, axis['key']), axis)
    row.update({
        'seed': params.get('SEED'),
        'cb_on': metadata.get('cb_on', cb_on_override),
        'low_res': metadata.get('low_res'),
        'num_cells': num_cells,
        'sim_duration_s': duration,
    })
    row.update(syn.summarize_run(results, num_cells, duration, t_start))
    row['source_path'] = os.path.abspath(path)
    return row


def build_table(run_dir: str, config: dict, cb_on_override: Optional[bool]) -> pd.DataFrame:
    """Rows for every raw file in run_dir. Files that fail to load are reported and skipped."""
    paths = find_raw_files(run_dir)
    if not paths:
        sys.exit(f"No jobs/*/output.pkl files under {run_dir}")

    rows, failed = [], []
    for n, path in enumerate(paths, start=1):
        print(f"[{n}/{len(paths)}] {path}")
        try:
            rows.append(job_row(path, config, cb_on_override))
        except Exception as err:  # keep going; report at the end
            failed.append((path, f'{type(err).__name__}: {err}'))
    for path, reason in failed:
        print(f"FAILED {path}: {reason}")
    if not rows:
        sys.exit("Every file failed; no table written")
    return pd.DataFrame(rows).sort_values(axis_keys(config) + ['seed'], ignore_index=True)


def report(table: pd.DataFrame, config: dict) -> None:
    """Print what each axis looks like, so a wrong config or a broken sweep is obvious."""
    print(f"\n{len(table)} rows")
    for key in axis_keys(config):
        values = sorted(float(v) for v in table[key].dropna().unique())
        print(f"  {key}: {len(values)} values {values}"
              + (f"  ({table[key].isna().sum()} rows NaN)" if table[key].isna().any() else ''))
    counts = table.groupby(axis_keys(config), dropna=False).size()
    print(f"  jobs per grid point: min {counts.min()}, max {counts.max()}")


def parse_bool(text: str) -> bool:
    if text.lower() in ('true', '1', 'yes'):
        return True
    if text.lower() in ('false', '0', 'no'):
        return False
    raise argparse.ArgumentTypeError(f"expected true/false, got {text!r}")


def main(argv=None):
    parser = argparse.ArgumentParser(description='Build results.csv from a sweep\'s raw output.')
    parser.add_argument('run_dir', help='Folder containing jobs/<job_id>/output.pkl')
    parser.add_argument('--config', help='sweep_config.yaml (default: <run_dir>/sweep_config.yaml)')
    parser.add_argument('--out', help='Output CSV (default: <run_dir>/results.csv)')
    parser.add_argument('--cb-on', type=parse_bool, default=None,
                        help='CB setting for old raw files that don\'t record it')
    args = parser.parse_args(argv)

    config_path = args.config or os.path.join(args.run_dir, 'sweep_config.yaml')
    try:
        config = load_sweep_config(config_path)
    except (OSError, ValueError) as err:
        sys.exit(f"Can't read sweep config {config_path}: {err}")

    table = build_table(args.run_dir, config, args.cb_on)
    out = args.out or os.path.join(args.run_dir, 'results.csv')
    table.to_csv(out, index=False)
    report(table, config)
    print(f"Wrote {out}")


if __name__ == '__main__':
    main()
