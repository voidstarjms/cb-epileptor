"""Generate one YAML per Condor job and a params_list.txt for Condor to queue over.

Each YAML is a complete, self-contained copy of the base params.yaml with the
sweep dimensions (CE, X0, Gintra, Ginter) overridden and a SEED field added.
Run once before submitting:

    python generate_params.py
"""
import os
import shutil
import argparse
import sys
from itertools import product

import numpy as np
import yaml

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_YAML = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..', 'parameters', 'params.yaml'))

def _write_param_value(out, param_name, param_val):
    #parameter = out[param_name]
    #if parameter is str:
    #parameter = out[param_name].split('*')
    #out[param_name] = f"{float(param_val)} *{parameter[1]}"
    #else:
        #out[param_name] = float(param_val)
    out[param_name] = param_val
    return

def _job_yaml(base, param1_val, param2_val, param3_val, param4_val, seed,
              param1_name = 'COUPLING_STRENGTH', param2_name = 'EPOP_BASE_EXCITE',
              param3_name = 'G_INTRA', param4_name = 'G_INTER'):
    """Return a base-copy dict with sweep dims and SEED overridden.

    Args:
        base (dict): Base YAML loaded from params.yaml.
        
        seed (int): Random seed for the job (stored as SEED field).

    Returns:
        dict: A copy of base with the four sweep dims and SEED replaced.
    """
    out = dict(base)
    # Override the constants that the model eqs actually read.
    _write_param_value(out, param1_name, param1_val)
    _write_param_value(out, param2_name, param2_val)
    _write_param_value(out, param3_name, param3_val)
    _write_param_value(out, param4_name, param4_val)
    # TODO Timed arrays are no longer updated here, so switching between use of the
    # timed arrays and the single values in the yaml needs to be handled in run.py
    # and model.py
    out['SEED'] = int(seed)
    return out

# TODO maybe remove this if no use is found
def _write_out_yaml(base, job, p1, p2, p3, p4, r, rel_paths):
    pass
    job_yaml = _job_yaml(base, p2, p1, p3, p4, seed=r)
    name = f'param_{job}.yaml'
    with open(os.path.join(OUT_DIR, name), 'w') as f:
        yaml.safe_dump(job_yaml, f, sort_keys=False)
    # params_list.txt holds paths relative to the Condor Initialdir
    # (which is set to src/sweep/ in setup_condor.sh).
    rel_paths.append(os.path.join('params', name))

def main():
    """Wipe params/, write one YAML per job in the Cartesian sweep grid.

    Also writes params_list.txt with one relative path per job, suitable for
    Condor's `Queue params_file from ...` directive (resolved against the
    Initialdir set in setup_condor.sh).
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', type=str, default='sweep_manifest.txt',
                        help='The path to the manifest file for the sweep')
    parser.add_argument('--granularity', '-g', type=int, default=4, required=True,
                        help='Number of steps along each parameter axis')
    parser.add_argument('--realizations', '-r', type=int, default=1,
                        help='Number of realizations for each parameter set')
    parser.add_argument('--outdir', type=str, default=os.path.join(SCRIPT_DIR,
                        'params', 'default'),
                        help='The directory to output the yaml files to')
    args = parser.parse_args()
    infile = args.infile
    n_realizations = args.realizations
    granularity = args.granularity
    WRAPPER_OUT_DIR = args.outdir
    OUT_DIR = os.path.join(WRAPPER_OUT_DIR, 'params')
    LIST_FILE = os.path.join(WRAPPER_OUT_DIR, 'params_list.txt')
    SIDECAR_FILE = os.path.join(WRAPPER_OUT_DIR, 'param_name_sidecar.txt')

    with open(BASE_YAML) as f:
        base = yaml.safe_load(f)

    # Read the sweep manifest
    params = []
    param_units = []
    param_mins = []
    param_maxs = []
    param_names = []
    with open(infile, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break
            # Split each line into [param_name, sweep_min, sweep_max, param_title]
            line_elements = line.split(maxsplit=4)
            line_element_count = len(line_elements)

            # Check split results to ensure correct structure
            if line_element_count != 5:
                print(f"Sweep manifest {f.name} contains malformed entries")
                sys.exit(1)

            # Check parameter name in base yaml
            if line_elements[0] not in base:
                print(f"{line_elements[0]} in sweep manifest is not a parameter in the base yaml")
                sys.exit(2)
            
            params.append(line_elements[0])
            param_units.append(line_elements[1])
            param_mins.append(line_elements[2])
            param_maxs.append(line_elements[3])
            param_names.append(line_elements[4])
    
    # TODO Ensure either 2 or 4 parameters were read (only 4 supported right now)
    param_count = len(params)
    if param_count != 2 and param_count != 4:
        print(f"Manifest must provide either 2 or 4 parameters, but got {param_count}")
        sys.exit(3)

    # Create files and directories for experiment
    if not os.path.isdir('params'):
        os.makedirs('params')
    if not os.path.isdir(WRAPPER_OUT_DIR):
        os.makedirs(WRAPPER_OUT_DIR)
    param_yaml_dir = os.path.join(WRAPPER_OUT_DIR, 'params')
    if not os.path.isdir(param_yaml_dir):
        os.makedirs(param_yaml_dir)
    
    # Write parameter names to sidecar text file

    with open(SIDECAR_FILE, "w") as f:
        for l in param_names:
            f.writelines(l)

    # Generate parameter ranges
    param1_values = np.linspace(float(param_mins[0]), float(param_maxs[0]), granularity)
    param2_values = np.linspace(float(param_mins[1]), float(param_maxs[1]), granularity)
    if param_count == 4:
        param3_values = np.linspace(float(param_mins[2]), float(param_maxs[2]), granularity)
        param4_values = np.linspace(float(param_mins[3]), float(param_maxs[3]), granularity)

    # Wipe yaml output directory
    if os.path.isdir(OUT_DIR):
        shutil.rmtree(OUT_DIR)
    os.makedirs(OUT_DIR)

    # TODO This only supports 4 parameters
    rel_paths = []
    combos = product(param1_values, param2_values, param3_values, param4_values,
                     range(1, n_realizations + 1))
    for job, (p1, p2, p3, p4, r) in enumerate(combos, start=1):
        job_yaml = _job_yaml(base, str(p2)+"*"+param_units[1],
                             str(p1)+"*"+param_units[0],
                             str(p3)+"*"+param_units[2],
                             str(p4)+"*"+param_units[3], r, 
                             param1_name=params[1],
                             param2_name=params[0],
                             param3_name=params[2],
                             param4_name=params[3])
        name = f'param_{job}.yaml'
        with open(os.path.join(OUT_DIR, name), 'w') as f:
            yaml.safe_dump(job_yaml, f, sort_keys=False)
        # params_list.txt holds paths relative to the Condor Initialdir
        # (which is set to src/sweep/ in setup_condor.sh).
        rel_paths.append(os.path.join('params', param_yaml_dir, name))

    with open(LIST_FILE, 'w') as f:
        f.write('\n'.join(rel_paths) + '\n')

    print(f"Wrote {len(rel_paths)} param YAMLs to {OUT_DIR}")
    print(f"Wrote {LIST_FILE} to {WRAPPER_OUT_DIR}")
    print(f"Grid: {param_names[1]} x {param_names[0]} x {param_names[2]} x {param_names[3]} x realizations = "
          f"{len(param1_values)} x {len(param2_values)} x {len(param3_values)} x "
          f"{len(param4_values)} x {n_realizations}")


if __name__ == '__main__':
    main()
