"""Load a flat YAML params file."""
import yaml
from brian2 import *


def load_params(path):
    """Load a flat YAML and resolve brian2 expressions.

    Resolves expressions such as '1*msecond' or '[1,2]*uS' by eval'ing each
    value against the brian2 namespace.

    Args:
        path (str): Path to the YAML file.

    Returns:
        dict: Mapping each YAML key to its resolved value.
    """
    with open(path, 'r') as f:
        raw = yaml.load(f, yaml.UnsafeLoader)
    ns = {}
    for k, v in raw.items():
        try:
            ns[k] = eval(str(v))
        except NameError:
            ns[k] = v
    return ns
