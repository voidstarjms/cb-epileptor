"""Load a sweep_config.yaml: the file that says what a sweep varies.

    base_params: ../../parameters/params.yaml   # relative to this file
    realizations: 2
    cb_on: false
    low_res: true
    sweep_params:                                 # order = plot order: x, y, facet row, facet col
      - {key: G_INTRA, unit: uS, min: 0.5, max: 2.0, steps: 4, label: "$G_{intra}$ (µS)"}
      - ...

The config is returned as a plain dict, the same shape as the YAML, with
defaults filled in and each swept parameter given a 'values' array.
"""
import os

import numpy as np
import yaml

MAX_SWEEP_PARAMS = 4


def load_sweep_config(path: str) -> dict:
    """Read a sweep config, fill in defaults, and build each swept parameter's values."""
    with open(path) as f:
        config = yaml.safe_load(f)

    # set default values for simulation settings
    config.setdefault('realizations', 1)
    config.setdefault('cb_on', False)
    config.setdefault('low_res', True)

    if len(config['sweep_params']) > MAX_SWEEP_PARAMS:
        raise ValueError(f"{path}: at most {MAX_SWEEP_PARAMS} sweep_params, got {len(config['sweep_params'])}")

    for sweep_param in config['sweep_params']:
        # construct the array of values a sweep param will take
        sweep_param['values'] = np.linspace(sweep_param['min'], sweep_param['max'], sweep_param['steps'])

    # base_params is written relative to the config file.
    config['base_params'] = os.path.join(os.path.dirname(path), config['base_params'])

    return config


def validate_sweep_params(config: dict) -> None:
    """Every swept key must exist in the base params, or it would be added as a new, unused key.

    Only generate_params calls this: copies of the config saved with a run's
    output sit in a different folder, so their base_params path no longer resolves.
    """
    with open(config['base_params']) as f:
        base = yaml.safe_load(f)
    missing = [p['key'] for p in config['sweep_params'] if p['key'] not in base]
    if missing:
        raise ValueError(f"not keys in {config['base_params']}: {missing}")


def sweep_param_keys(config: dict) -> list:
    """The swept parameter keys, in plot order."""
    return [p['key'] for p in config['sweep_params']]
