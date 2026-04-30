"""Aggregation script — run after all Condor jobs finish.

Reads all per-job result pkls under data/results/, groups by (ce, x0, Gintra,
Ginter), computes mean and SD chi across realizations, and writes heatmap PNGs
to figures/. The run number is prefixed onto output filenames so successive
sweeps don't overwrite each other.

Usage:

    python aggregate.py             # default 2-D heatmap
    python aggregate.py --full      # multifaceted (Gintra x Ginter) panel grid
"""
import os
import sys
import argparse
import pickle
import numpy as np
from collections import defaultdict

sys.path.append("..")
from run import FilePaths
import plotting.sync_heatmap as ps

results_dir = os.path.join('data', 'results')

if not os.path.exists(results_dir):
    print(f"No results directory found at {results_dir}")
    sys.exit(1)

# Load all job result pkls
all_results = []
for fname in sorted(os.listdir(results_dir)):
    if fname.endswith('.pkl'):
        with open(os.path.join(results_dir, fname), 'rb') as f:
            all_results.append(pickle.load(f))

if len(all_results) == 0:
    print("No result files found in data/results/")
    sys.exit(1)

print(f"Loaded {len(all_results)} job results")

parser = argparse.ArgumentParser(description='Aggregate chi values from experiments into a heat map')
parser.add_argument('--full', default=False, action='store_true', help='Plot multifaceted heatmap to include synaptic conductances')
args = parser.parse_args()

# Reconstruct param axes from results
ce_values = sorted(set(round(r['ce'], 6) for r in all_results))
x0_values = sorted(set(round(r['x0'], 6) for r in all_results))
Gintra_values = sorted(set(round(r['Gintra'], 6) for r in all_results))
Ginter_values = sorted(set(round(r['Ginter'], 6) for r in all_results))
print(f"CE values  ({len(ce_values)}): {[round(v,3) for v in ce_values]}")
print(f"X0 values ({len(x0_values)}): {[round(v,3) for v in x0_values]}")
print(f"Gintra values ({len(Gintra_values)}): {[round(v,3) for v in Gintra_values]}")
print(f"Ginter values ({len(Ginter_values)}): {[round(v,3) for v in Ginter_values]}")

p1_label = r'$C_E$ (coupling strength)'
p2_label = r'$x_0$ (epileptogenicity)'

# Group chi values by (ce, x0, Gintra, Ginter) across realizations
chi_by_point = defaultdict(list)
for r in all_results:
    key = (round(r['ce'], 6), round(r['x0'], 6), round(r['Gintra'], 6), round(r['Ginter'], 6))
    chi_by_point[key].append(r['chi'])

figures_dir = 'figures'
os.makedirs(figures_dir, exist_ok=True)
run_num_file = 'current_run.txt'
filepaths = FilePaths(
    None, 'figures'
)

# 2-axis plot
if args.full == False:
    # Build mean and SD grids — shape: (len(x0), len(ce))
    chi_grid = np.full((len(x0_values), len(ce_values)), np.nan)
    chi_sd   = np.full((len(x0_values), len(ce_values)), np.nan)
    Gintra_min = min(Gintra_values)
    Ginter_min = min(Ginter_values)
    for (ce, x0, Gintra, Ginter), chis in chi_by_point.items():
        # Only generate chi_grid elements for min G values for the simple graph
        if Gintra != Gintra_min or Ginter != Ginter_min:
            continue
        i = ce_values.index(ce)
        j = x0_values.index(x0)
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
                  np.array(ce_values), np.array(x0_values),
                  p1_label, p2_label,
                  run_num=run_num)

    print(f"Saved: {run_num}_synchrony_chi_mean.png and {run_num}_synchrony_chi_sd.png")
else:
    p3_label = r'$G_{intra}$ (intrapop conductance)'
    p4_label = r'$G_{inter}$ (interpop conductance)'

    # Build mean and SD grids — shape: (len(x0), len(ce), len(Gintra), len(Ginter))
    chi_grid = np.full((len(x0_values), len(ce_values), len(Gintra_values), len(Ginter_values)), np.nan)
    chi_sd   = np.full((len(x0_values), len(ce_values), len(Gintra_values), len(Ginter_values)), np.nan)
    for (ce, x0, Gintra, Ginter), chis in chi_by_point.items():
        i = ce_values.index(ce)
        j = x0_values.index(x0)
        k = Gintra_values.index(Gintra)
        l = Ginter_values.index(Ginter)
        chi_grid[j, i, k, l] = np.mean(chis)
        chi_sd[j, i, k, l]   = np.std(chis)

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
        while os.path.exists(os.path.join(figures_dir, f'{run_num}_sweep_debug_full')):
            run_num += 1

    # TODO This doesn't compute STD yet
    ps.plot_synchrony_multifaceted(filepaths, chi_grid,
                  np.array(ce_values), np.array(x0_values), np.array(Gintra_values),
                  np.array(Ginter_values), p1_label, p2_label, p3_label, p4_label)
    
    print("Saved: {run_num}_full_param_space.png")
    #print(f"Saved: {run_num}_synchrony_chi_mean.png and {run_num}_synchrony_chi_sd.png")
