"""Load a flat YAML params file — mirrors Grant's brian_utils/model_IO.py."""
import yaml
from brian2 import *


def load_params(path):
    with open(path, 'r') as f:
        raw = yaml.load(f, yaml.UnsafeLoader)
    ns = {}
    for k, v in raw.items():
        try:
            ns[k] = eval(str(v))
        except NameError:
            ns[k] = v
    return ns
