"""Aggregation script — run after all Condor jobs finish.

Reads all per-job result pkls under data/results/, groups by (ce, x0, Gintra,
Ginter), computes mean and SD chi across realizations, and writes heatmap PNGs
to figures/. The run number is prefixed onto output filenames so successive
sweeps don't overwrite each other.

Usage:

    python aggregate.py             # default 2D heatmap
    python aggregate.py --full      # multifaceted (4D) panel grid
"""
import os
import sys
import argparse
import pickle
import numpy as np
from collections import defaultdict
from brian2 import Quantity

sys.path.append("..")
from run import FilePaths
import plotting.sync_heatmap as ps

parser = argparse.ArgumentParser(description='Aggregate synchrony values from experiments into a heat map')
parser.add_argument('--full', default=False, action='store_true', help='Plot multifaceted heatmap to include synaptic conductances')
parser.add_argument('--name', default='default', type=str, help='Name of run to load. Default: default')
args = parser.parse_args()
run_name = args.name

param_path = os.path.join('params', run_name)
data_path = os.path.join('..', '..', 'output', 'sweep', run_name)
results_dir = os.path.join(data_path, 'data', 'results')
sidecar_file = os.path.join(param_path, 'param_name_sidecar.txt')

if not os.path.exists(param_path):
    print(f"Run parameter directory with name {run_name} not found")
    sys.exit(1)

if not os.path.exists(data_path):
    print(f"Run output directory with name {run_name} not found")
    sys.exit(1)

if not os.path.exists(results_dir):
    print(f"No results directory found at {results_dir}")
    sys.exit(1)

if not os.path.exists(sidecar_file):
    print(f"Run {run_name} is missing parameter sidecar file")
    sys.exit(1)

# Open sidecar file for run to get parameters and their display names
params = []
param_names = []
with open(sidecar_file, "r") as f:
    while True:
        line = f.readline()
        if line == "":
            break
        elements = line.split(maxsplit=1)
        params.append(elements[0])
        param_names.append(elements[1].strip())

# Load all job result pkls
all_results = []
for fname in sorted(os.listdir(results_dir)):
    if fname.endswith('.pkl'):
        with open(os.path.join(results_dir, fname), 'rb') as f:
            contents = pickle.load(f)
            all_results.append(contents)

if len(all_results) == 0:
    print(f"No result files found in {results_dir}")
    sys.exit(1)

print(f"Loaded {len(all_results)} job results")

def _strip_units(x):
    if type(x) == Quantity:
        u = x.get_best_unit()
        return x / u
    else:
        return x

# Reconstruct param axes from results
# TODO This can definitely be done with an array
param1_vals = sorted(set(round(_strip_units(r['params'][params[0]]), 6) for r in all_results))
param2_vals = sorted(set(round(_strip_units(r['params'][params[1]]), 6) for r in all_results))
param3_vals = sorted(set(round(_strip_units(r['params'][params[2]]), 6) for r in all_results))
param4_vals = sorted(set(round(_strip_units(r['params'][params[3]]), 6) for r in all_results))
print(f"{param_names[0]} values  ({len(param1_vals)}): {[round(v, 3) for v in param1_vals]}")
print(f"{param_names[1]} values ({len(param2_vals)}): {[round(v, 3) for v in param2_vals]}")
print(f"{param_names[2]} values ({len(param3_vals)}): {[round(v, 3) for v in param3_vals]}")
print(f"{param_names[3]} values ({len(param4_vals)}): {[round(v, 3) for v in param4_vals]}")

p1_label = param_names[0]
p2_label = param_names[1]

# Group chi values by (ce, x0, Gintra, Ginter) across realizations
chi_by_point = defaultdict(list)
for r in all_results:
    #print(r)
    key = (round(_strip_units(r['params'][params[0]]), 6),
            round(_strip_units(r['params'][params[1]]), 6),
            round(_strip_units(r['params'][params[2]]), 6),
            round(_strip_units(r['params'][params[3]]), 6))
    chi_by_point[key].append(r['r_exc'] + r['r_inh'])

figures_dir = os.path.join('..', '..', 'figures', 'sweep')
os.makedirs(figures_dir, exist_ok=True)
run_num_file = 'current_run.txt'
filepaths = FilePaths(
    None, 'figures'
)

# 2-axis plot
if args.full == False:
    # Build mean and SD grids — shape: (len(x0), len(ce))
    chi_grid = np.full((len(param2_vals), len(param1_vals)), np.nan)
    chi_sd   = np.full((len(param2_vals), len(param1_vals)), np.nan)
    p3_min = min(param3_vals)
    p4_min = min(param4_vals)
    for (p1, p2, p3, p4), chis in chi_by_point.items():
        # Only generate chi_grid elements for min G values for the simple graph
        if p3 != p3_min or p4 != p4_min:
            continue
        i = param2_vals.index(p2)
        j = param1_vals.index(p1)
        chi_grid[j, i] = np.mean(chis)
        chi_sd[j, i]   = np.std(chis)

    missing = int(np.sum(np.isnan(chi_grid)))
    if missing > 0:
        print(f"Warning: {missing} grid points are missing results")
    else:
        print("All grid points complete.")

    # Read run number set by setup_condor.sh so logs/figures/debug all share the same number
    if os.path.exists(run_num_file):
        with open(run_num_file) as f:
            run_num = int(f.read().strip())
    else:
        # Fallback: pick next available number
        run_num = 1
        while os.path.exists(os.path.join(figures_dir, f'{run_num}_sweep_debug')):
            run_num += 1


    ps.plot_synchrony(filepaths, chi_grid, chi_sd,
                  np.array(param1_vals), np.array(param2_vals),
                  p1_label, p2_label, run_name)

# 4-axis plot
else:
    p3_label = param_names[2]
    p4_label = param_names[3]

    # Build mean and SD grids — shape: (len(x0), len(ce), len(Gintra), len(Ginter))
    chi_grid = np.full((len(param2_vals), len(param1_vals), len(param3_vals), len(param4_vals)), np.nan)
    chi_sd   = np.full((len(param2_vals), len(param1_vals), len(param3_vals), len(param4_vals)), np.nan)
    for (p1, p2, p3, p4), chis in chi_by_point.items():
        i = param2_vals.index(p2)
        j = param1_vals.index(p1)
        l = param3_vals.index(p3)
        k = param4_vals.index(p4)
        chi_grid[i, j, k, l] = np.mean(chis)
        chi_sd[i, j, k, l]   = np.std(chis)
        print(p2, p1, p3, p4, chi_grid[i, j, k, l])

    missing = int(np.sum(np.isnan(chi_grid)))
    if missing > 0:
        print(f"Warning: {missing} grid points are missing results")
    else:
        print("All grid points complete.")

    # Read run number set by setup_condor.sh so logs/figures/debug all share the same number
   
    #if os.path.exists(run_num_file):
    #    with open(run_num_file) as f:
    #        run_num = int(f.read().strip())
    #else:
    #    # Fallback: pick next available number
    #    run_num = 1
    #    while os.path.exists(os.path.join(figures_dir, f'{run_num}_sweep_debug_full')):
    #        run_num += 1

    ps.plot_synchrony_multifaceted(filepaths, chi_grid, chi_sd,
                  np.array(param2_vals), np.array(param1_vals), np.array(param4_vals),
                  np.array(param3_vals), p2_label, p1_label, p3_label, p4_label,
                  run_name)
