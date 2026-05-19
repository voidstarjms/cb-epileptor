"""Generate one YAML per Condor job and a params_list.txt for Condor to queue over.

Each YAML is a complete, self-contained copy of the base params.yaml with the
sweep dimensions (CE, X0, Gintra, Ginter) overridden and a SEED field added.
Run once before submitting:

    python generate_params.py
"""
import os
import shutil
from itertools import product

import numpy as np
import yaml

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_YAML = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..', 'params.yaml'))
OUT_DIR = os.path.join(SCRIPT_DIR, 'params')
LIST_FILE = os.path.join(SCRIPT_DIR, 'params_list.txt')

GRANULARITY = 2
ce_values     = np.linspace(0.0, 0.1, GRANULARITY)    # CE (coupling strength)
x0_values     = np.linspace(-4.5, -2.5, GRANULARITY)  # X0 (epileptogenicity)
gintra_values = np.linspace(0.5, 2.0, GRANULARITY)    # Gintra (intra-pop conductance, uS)
ginter_values = np.linspace(1.0, 4.0, GRANULARITY)    # Ginter (inter-pop conductance, uS)
n_realizations = 5


def _job_yaml(base, ce, x0, gintra, ginter, seed):
    """Return a base-copy dict with sweep dims and SEED overridden.

    Args:
        base (dict): Base YAML loaded from params.yaml.
        ce (float): Coupling strength.
        x0 (float): Epileptogenicity (HR_X_NAUGHT).
        gintra (float): Intra-pop conductance in microsiemens.
        ginter (float): Inter-pop conductance in microsiemens.
        seed (int): Random seed for the job (stored as SEED field).

    Returns:
        dict: A copy of base with the four sweep dims and SEED replaced.
    """
    out = dict(base)
    # Override the constants that the model eqs actually read.
    out['COUPLING_STRENGTH'] = float(ce)
    out['HR_X_NAUGHT']       = float(x0)
    out['G_INTRA']           = f'{float(gintra)} * uS'
    out['G_INTER']           = f'{float(ginter)} * uS'
    # Also override the timed-array schedules so plots/labels stay consistent.
    out['COUPLING_VALS']     = [float(ce)]
    out['X_NAUGHT_VALS']     = [float(x0)]
    out['G_INTRA_VALS']      = f'[{float(gintra)}] * uS'
    out['G_INTER_VALS']      = f'[{float(ginter)}] * uS'
    out['SEED']              = int(seed)
    return out


def main():
    """Wipe params/, write one YAML per job in the Cartesian sweep grid.

    Also writes params_list.txt with one relative path per job, suitable for
    Condor's `Queue params_file from ...` directive (resolved against the
    Initialdir set in setup_condor.sh).
    """
    with open(BASE_YAML) as f:
        base = yaml.safe_load(f)

    if os.path.isdir(OUT_DIR):
        shutil.rmtree(OUT_DIR)
    os.makedirs(OUT_DIR)

    rel_paths = []
    combos = product(x0_values, ce_values, gintra_values, ginter_values,
                     range(1, n_realizations + 1))
    for job, (x0, ce, gintra, ginter, r) in enumerate(combos, start=1):
        job_yaml = _job_yaml(base, ce, x0, gintra, ginter, seed=r)
        name = f'param_{job}.yaml'
        with open(os.path.join(OUT_DIR, name), 'w') as f:
            yaml.safe_dump(job_yaml, f, sort_keys=False)
        # params_list.txt holds paths relative to the Condor Initialdir
        # (which is set to src/sweep/ in setup_condor.sh).
        rel_paths.append(os.path.join('params', name))

    with open(LIST_FILE, 'w') as f:
        f.write('\n'.join(rel_paths) + '\n')

    print(f"Wrote {len(rel_paths)} param YAMLs to {OUT_DIR}")
    print(f"Wrote {LIST_FILE}")
    print(f"Grid: CE x X0 x Gintra x Ginter x realizations = "
          f"{len(ce_values)} x {len(x0_values)} x {len(gintra_values)} x "
          f"{len(ginter_values)} x {n_realizations}")


if __name__ == '__main__':
    main()
