"""Load a sweep_config.yaml: the file that says what a sweep varies.

    base_params: ../../parameters/params.yaml   # relative to this file
    realizations: 2
    cb_on: false
    low_res: true
    axes:                                         # order = plot order: x, y, facet row, facet col
      - {key: G_INTRA, unit: uS, min: 0.5, max: 2.0, steps: 4, label: "$G_{intra}$ (µS)"}
      - ...

The config is returned as a plain dict, the same shape as the YAML, with
defaults filled in and each axis given a 'values' array.

generate_params.py needs the full file. build_table.py and plot_sweep.py only
need each axis's key, unit and label, so a hand-written config for an old sweep
can leave out min/max/steps and base_params.
"""
import os

import numpy as np
import yaml

MAX_AXES = 4


def load_sweep_config(path: str) -> dict:
    """Read a sweep config, fill in defaults, and build each axis's values."""
    with open(path) as f:
        config = yaml.safe_load(f)

    config.setdefault('realizations', 1)
    config.setdefault('cb_on', False)
    config.setdefault('low_res', True)

    if len(config['axes']) > MAX_AXES:
        raise ValueError(f"{path}: at most {MAX_AXES} axes, got {len(config['axes'])}")

    for axis in config['axes']:
        if 'min' in axis:  # configs for old sweeps have no range
            axis['values'] = np.linspace(axis['min'], axis['max'], axis['steps'])

    if 'base_params' in config:
        # base_params is written relative to the config file.
        config['base_params'] = os.path.join(os.path.dirname(path), config['base_params'])

    return config


def check_keys_in_base(config: dict) -> None:
    """Every swept key must exist in the base params, or it would be added as a new, unused key.

    Only generate_params calls this: copies of the config saved with a run's
    output sit in a different folder, so their base_params path no longer resolves.
    """
    with open(config['base_params']) as f:
        base = yaml.safe_load(f)
    missing = [axis['key'] for axis in config['axes'] if axis['key'] not in base]
    if missing:
        raise ValueError(f"not keys in {config['base_params']}: {missing}")


def axis_keys(config: dict) -> list:
    """The swept parameter keys, in plot order."""
    return [axis['key'] for axis in config['axes']]
